#include <cmath>
#include "Uno.hpp"
#include "Iterate.hpp"
#include "Logger.hpp"
#include "Statistics.hpp"
#include "Timer.hpp"
#include "Preprocessing.hpp"

Uno::Uno(GlobalizationMechanism& globalization_mechanism, double tolerance, int max_iterations) : globalization_mechanism(
      globalization_mechanism), tolerance(tolerance), max_iterations(max_iterations) {
}

Result Uno::solve(Problem& problem, std::vector<double>& x, Multipliers& multipliers, bool preprocessing) {
   Timer timer{};
   timer.start();
   int major_iterations = 0;

   INFO << "Problem " << problem.name << "\n";
   INFO << problem.number_variables << " variables, " << problem.number_constraints << " constraints\n";
   INFO << "Problem type: " << Problem::type_to_string[problem.type] << "\n";

   /* project x into the bounds */
   Subproblem::project_point_in_bounds(x, problem.variables_bounds);
   if (preprocessing) {
      /* preprocessing phase: satisfy linear constraints */
      Preprocessing::apply(problem, x, multipliers);
   }

   Statistics statistics = Uno::create_statistics();
   /* use the current point to initialize the strategies and generate the initial point */
   Iterate current_iterate = this->globalization_mechanism.initialize(statistics, problem, x, multipliers);

   TerminationStatus termination_status = NOT_OPTIMAL;
   try {
      /* check for convergence */
      while (!this->termination_criterion_(termination_status, major_iterations)) {
         statistics.new_line();
         major_iterations++;
         DEBUG << "Current iterate\n" << current_iterate << "\n";

         /* compute an acceptable iterate by solving a subproblem at the current point */
         auto [new_iterate, direction] = this->globalization_mechanism.compute_acceptable_iterate(statistics, problem, current_iterate);
         DEBUG << "Next iterate\n" << new_iterate;

         Uno::add_statistics(statistics, new_iterate, major_iterations);
         statistics.print_current_line();

         // compute the status of the new iterate
         termination_status = this->check_termination(problem, new_iterate, direction.norm, direction.objective_multiplier);
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
   timer.stop();

   int number_subproblems_solved = this->globalization_mechanism.constraint_relaxation_strategy.subproblem.number_subproblems_solved;
   Result result =
         {termination_status, current_iterate, problem.number_variables, problem.number_constraints, major_iterations, timer.get_time(),
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

void Uno::add_statistics(Statistics& statistics, const Iterate& new_iterate, int major_iterations) {
   statistics.add_statistic("major", major_iterations);
   statistics.add_statistic("f", new_iterate.objective);
   statistics.add_statistic("||c||", new_iterate.residuals.constraints);
   statistics.add_statistic("complementarity", new_iterate.residuals.complementarity);
   statistics.add_statistic("KKT", new_iterate.residuals.KKT);
   statistics.add_statistic("FJ", new_iterate.residuals.FJ);
}

bool Uno::termination_criterion_(TerminationStatus current_status, int iteration) const {
   return current_status != NOT_OPTIMAL || this->max_iterations <= iteration;
}

TerminationStatus
Uno::check_termination(Problem& problem, Iterate& current_iterate, double step_norm, double objective_multiplier) const {
   TerminationStatus status = NOT_OPTIMAL;

   if (current_iterate.residuals.complementarity <= this->tolerance * (double) (current_iterate.x.size() + problem.number_constraints)) {
      // feasible and KKT point
      if (current_iterate.residuals.KKT <= this->tolerance * std::sqrt(current_iterate.x.size())) {
         if (current_iterate.residuals.constraints <= this->tolerance * (double) current_iterate.x.size()) {
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

   INFO << "Feasibility measure:\t\t" << this->solution.progress.feasibility << "\n";
   INFO << "Optimality measure:\t\t" << this->solution.progress.objective << "\n";

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
