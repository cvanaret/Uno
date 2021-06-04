#include <iostream>
#include <ctime>
#include <cmath>
#include "Uno.hpp"
#include "Iterate.hpp"
#include "Logger.hpp"
#include "BQPDSolver.hpp"
#include "Statistics.hpp"

Uno::Uno(GlobalizationMechanism& globalization_mechanism, double tolerance, int max_iterations) :
globalization_mechanism(globalization_mechanism), tolerance(tolerance), max_iterations(max_iterations) {
}

Result Uno::solve(Problem& problem, std::vector<double>& x, Multipliers& multipliers, bool preprocessing) {
   std::clock_t c_start = std::clock();
   int major_iterations = 0;

   INFO << "Problem " << problem.name << "\n";
   INFO << problem.number_variables << " variables, " << problem.number_constraints << " constraints\n";
   INFO << "Problem type: " << Problem::type_to_string[problem.type] << "\n";

   /* project x into the bounds */
   Subproblem::project_point_in_bounds(x, problem.variables_bounds);

   if (preprocessing) {
      /* preprocessing phase: satisfy linear constraints */
      this->preprocessing(problem, x, multipliers);
   }

   Statistics statistics = Uno::create_statistics();
   /* use the current point to initialize the strategies and generate the initial point */
   Iterate current_iterate = this->globalization_mechanism.initialize(statistics, problem, x, multipliers);
   DEBUG << "Initial iterate\n" << current_iterate << "\n";

   TerminationStatus termination_status = NOT_OPTIMAL;
   try {
      /* check for convergence */
      while (!this->termination_criterion_(termination_status, major_iterations)) {
         statistics.new_line();

         major_iterations++;
         DEBUG << "\n\t\tUNO iteration " << major_iterations << "\n";
         DEBUG << "Current point: ";
         print_vector(DEBUG, current_iterate.x);
         /* update the current point */
         auto[new_iterate, direction] = this->globalization_mechanism.compute_acceptable_iterate(statistics, problem, current_iterate);

         statistics.add_statistic("major", major_iterations);
         statistics.add_statistic("f", new_iterate.objective);
         statistics.add_statistic("||c||", new_iterate.residuals.constraints);
         statistics.add_statistic("complementarity", new_iterate.residuals.complementarity);
         statistics.add_statistic("KKT", new_iterate.residuals.KKT);
         statistics.add_statistic("FJ", new_iterate.residuals.FJ);
         statistics.print_current_line();
         DEBUG << "Next iterate\n" << new_iterate;
         termination_status = this->compute_termination_status_(problem, new_iterate, direction.norm, direction.objective_multiplier);
         current_iterate = new_iterate;
      }
   }
   catch (std::invalid_argument& exception) {
      ERROR << exception.what();
   }
   catch (std::runtime_error& exception) {
      ERROR << exception.what();
   }
   statistics.print_footer();
   std::clock_t c_end = std::clock();
   double cpu_time = (c_end - c_start) / (double) CLOCKS_PER_SEC;

   int number_subproblems_solved = this->globalization_mechanism.globalization_strategy.subproblem.number_subproblems_solved;
   Result result =
         {termination_status, current_iterate, problem.number_variables, problem.number_constraints, major_iterations, cpu_time,
          Iterate::number_eval_objective, Iterate::number_eval_constraints, Iterate::number_eval_jacobian, Iterate::number_eval_hessian,
          number_subproblems_solved};
   return result;
}

Statistics Uno::create_statistics() {
   Statistics statistics;
   statistics.add_column("major", Statistics::int_width, 1);
   statistics.add_column("minor", Statistics::int_width, 2);
   statistics.add_column("step norm", Statistics::double_width, 31);
   statistics.add_column("f", Statistics::double_width, 100);
   statistics.add_column("||c||", Statistics::double_width, 101);
   statistics.add_column("complementarity", Statistics::double_width, 104);
   statistics.add_column("KKT", Statistics::double_width, 105);
   statistics.add_column("FJ", Statistics::double_width, 106);
   return statistics;
}

bool Uno::termination_criterion_(TerminationStatus current_status, int iteration) const {
   return current_status != NOT_OPTIMAL || this->max_iterations <= iteration;
}

TerminationStatus
Uno::compute_termination_status_(Problem& problem, Iterate& current_iterate, double step_norm, double objective_multiplier) const {
   TerminationStatus status = NOT_OPTIMAL;

   if (current_iterate.residuals.complementarity <= this->tolerance * (double) (current_iterate.x.size() + problem.number_constraints)) {
      // feasible and KKT point
      if (current_iterate.residuals.constraints <= this->tolerance * (double) current_iterate.x.size()) {
         if (current_iterate.residuals.KKT <= this->tolerance * std::sqrt(current_iterate.x.size())) {
            status = KKT_POINT;
         }
      }
         // infeasible and FJ point
      else if (current_iterate.residuals.FJ <= this->tolerance * std::sqrt(current_iterate.x.size())) {
         status = FJ_POINT;
      }
   }
   else if (step_norm <= this->tolerance / 100.) {
      if (current_iterate.residuals.constraints <= this->tolerance * (double) current_iterate.x.size()) {
         status = FEASIBLE_SMALL_STEP;
      }
      else {
         status = INFEASIBLE_SMALL_STEP;
      }
   }

   // if convergence, correct the multipliers
   if (status != NOT_OPTIMAL && 0. < objective_multiplier) {
      for (double& multiplier_j: current_iterate.multipliers.constraints) {
         multiplier_j /= objective_multiplier;
      }
      for (size_t i = 0; i < current_iterate.x.size(); i++) {
         current_iterate.multipliers.lower_bounds[i] /= objective_multiplier;
         current_iterate.multipliers.upper_bounds[i] /= objective_multiplier;
      }
   }
   return status;
}

void Uno::preprocessing(Problem& problem, std::vector<double>& x, Multipliers& multipliers) {
   /* linear constraints */
   INFO << "Preprocessing phase: the problem has " << problem.linear_constraints.size() << " linear constraints\n";
   if (0 < problem.linear_constraints.size()) {
      std::vector<double> constraints = problem.evaluate_constraints(x);

      int infeasible_linear_constraints = 0;
      for (std::pair<int, int> element: problem.linear_constraints) {
         int j = element.first;
         if (constraints[j] < problem.constraint_bounds[j].lb || problem.constraint_bounds[j].ub < constraints[j]) {
            infeasible_linear_constraints++;
         }
      }
      INFO << "There are " << infeasible_linear_constraints << " infeasible linear constraints at the initial point\n";

      if (0 < infeasible_linear_constraints) {
         INFO << "Current point: ";
         print_vector(INFO, x);
         int number_constraints = (int) problem.linear_constraints.size();
         BQPDSolver solver(problem.number_variables, number_constraints, problem.number_variables, true);

         int fortran_indexing = 1;
         CSCMatrix hessian = CSCMatrix::identity(problem.number_variables, fortran_indexing);
         SparseGradient linear_objective; // empty
         std::vector<double> d0(problem.number_variables);
         // constraints Jacobian
         std::vector<SparseGradient> constraints_jacobian(number_constraints);
         for (const auto[j, linear_constraint_index]: problem.linear_constraints) {
            problem.constraint_gradient(x, j, constraints_jacobian[linear_constraint_index]);
         }
         // variables bounds
         std::vector<Range> variables_bounds(problem.number_variables);
         for (size_t i = 0; i < problem.number_variables; i++) {
            variables_bounds[i] = {problem.variables_bounds[i].lb - x[i], problem.variables_bounds[i].ub - x[i]};
         }
         // constraints bounds
         std::vector<Range> constraints_bounds(number_constraints);
         for (const auto[j, linear_constraint_index]: problem.linear_constraints) {
            constraints_bounds[linear_constraint_index] =
                  {problem.constraint_bounds[j].lb - constraints[j], problem.constraint_bounds[j].ub - constraints[j]};
         }
         Direction direction = solver.solve_QP(variables_bounds, constraints_bounds, linear_objective, constraints_jacobian, hessian, d0);
         if (direction.status == INFEASIBLE) {
            throw std::runtime_error("Linear constraints cannot be satisfied");
         }

         std::vector<double> feasible_x = add_vectors(x, direction.x, 1.);
         x = feasible_x;
         // copy bound multipliers
         multipliers.lower_bounds = direction.multipliers.lower_bounds;
         multipliers.upper_bounds = direction.multipliers.upper_bounds;
         // copy constraint multipliers
         for (const auto[j, linear_constraint_index]: problem.linear_constraints) {
            multipliers.constraints[j] = direction.multipliers.constraints[linear_constraint_index];
         }
         INFO << "Linear feasible initial point: ";
         print_vector(INFO, x);
         INFO << "\n";
      }
   }
}

void Result::display(bool print_solution) {
   INFO << "\n";
   INFO << "UNO v1: optimization summary\n";
   INFO << "==============================\n";

   INFO << "Status:\t\t\t\t";
   if (this->status == KKT_POINT) {
      INFO << "Converged with KKT point\n";
   }
   else if (this->status == FJ_POINT) {
      INFO << "Converged with FJ point\n";
   }
   else if (this->status == FEASIBLE_SMALL_STEP) {
      INFO << "Converged with feasible small step\n";
   }
   else if (this->status == INFEASIBLE_SMALL_STEP) {
      INFO << "Converged with infeasible small step\n";
   }
   else { // NOT_OPTIMAL
      INFO << "Irregular termination\n";
   }

   INFO << "Objective value:\t\t" << this->solution.objective << "\n";
   INFO << "Constraint residual:\t\t" << this->solution.residuals.constraints << "\n";
   INFO << "KKT residual:\t\t\t" << this->solution.residuals.KKT << "\n";
   INFO << "FJ residual:\t\t\t" << this->solution.residuals.FJ << "\n";
   INFO << "Complementarity residual:\t" << this->solution.residuals.complementarity << "\n";

   INFO << "Feasibility measure:\t\t" << this->solution.feasibility_measure << "\n";
   INFO << "Optimality measure:\t\t" << this->solution.optimality_measure << "\n";

   if (print_solution) {
      INFO << "Primal solution:\t\t";
      print_vector(INFO, this->solution.x);
      INFO << "Lower bound multipliers:\t";
      print_vector(INFO, this->solution.multipliers.lower_bounds);
      INFO << "Upper bound multipliers:\t";
      print_vector(INFO, this->solution.multipliers.upper_bounds);
      INFO << "Constraint multipliers:\t\t";
      print_vector(INFO, this->solution.multipliers.constraints);
   }

   INFO << "CPU time:\t\t\t" << this->cpu_time << "s\n";
   INFO << "Iterations:\t\t\t" << this->iteration << "\n";
   INFO << "Objective evaluations:\t\t" << this->objective_evaluations << "\n";
   INFO << "Constraints evaluations:\t" << this->constraint_evaluations << "\n";
   INFO << "Jacobian evaluations:\t\t" << this->jacobian_evaluations << "\n";
   INFO << "Hessian evaluations:\t\t" << this->hessian_evaluations << "\n";
   INFO << "Number of subproblems solved:\t" << this->number_subproblems_solved << "\n";
}
