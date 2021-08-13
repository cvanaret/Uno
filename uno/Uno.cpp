#include <cmath>
#include "Uno.hpp"
#include "optimization_problem/Iterate.hpp"
#include "optimization_problem/Preprocessing.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"
#include "tools/Timer.hpp"

Uno::Uno(GlobalizationMechanism& globalization_mechanism, double tolerance, int max_iterations) : globalization_mechanism(
      globalization_mechanism), tolerance(tolerance), max_iterations(max_iterations) {
}

Result Uno::solve(const Problem& problem, std::vector<double>& x, Multipliers& multipliers, bool use_preprocessing) {
   Timer timer{};
   timer.start();
   int major_iterations = 0;

   std::cout << "\nProblem " << problem.name << "\n";
   std::cout << problem.number_variables << " variables, " << problem.number_constraints << " constraints\n";
   std::cout << "Problem type: " << Problem::type_to_string[problem.type] << "\n";

   /* project x into the bounds */
   problem.project_point_in_bounds(x);
   if (use_preprocessing) {
      /* preprocessing phase: satisfy linear constraints */
      Preprocessing::apply(problem, x, multipliers);
   }

   Statistics statistics = Uno::create_statistics();
   /* use the current point to initialize the strategies and generate the initial iterate */
   Iterate current_iterate = this->globalization_mechanism.initialize(statistics, problem, x, multipliers);

   TerminationStatus termination_status = NOT_OPTIMAL;
   try {
      /* check for convergence */
      while (!this->termination_criterion(termination_status, major_iterations)) {
         statistics.new_line();
         major_iterations++;
         DEBUG << "\n########## Outer iteration " << major_iterations << "\n";
         DEBUG << "Current iterate\n" << current_iterate << "\n";

         /* compute an acceptable iterate by solving a subproblem at the current point */
         auto [new_iterate, direction_norm, objective_multiplier] = this->globalization_mechanism.compute_acceptable_iterate(statistics, problem, current_iterate);

         Uno::add_statistics(statistics, new_iterate, major_iterations);
         if (Logger::logger_level == INFO) statistics.print_current_line();

         // compute the status of the new iterate
         termination_status = this->check_termination(problem, new_iterate, direction_norm, objective_multiplier);
         current_iterate = std::move(new_iterate);
      }
   }
   catch (std::invalid_argument& exception) {
      ERROR << exception.what();
   }
   catch (std::runtime_error& exception) {
      ERROR << exception.what();
   }
   if (Logger::logger_level == INFO) statistics.print_footer();
   timer.stop();

   int number_subproblems_solved = this->globalization_mechanism.get_number_subproblems_solved();
   int hessian_evaluation_count = this->globalization_mechanism.get_hessian_evaluation_count();
   Result result =
         {termination_status, std::move(current_iterate), problem.number_variables, problem.number_constraints, major_iterations, timer.get_time(),
          Iterate::number_eval_objective, Iterate::number_eval_constraints, Iterate::number_eval_jacobian, hessian_evaluation_count,
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
   statistics.add_statistic("||c||", new_iterate.errors.constraints);
   statistics.add_statistic("complementarity", new_iterate.errors.complementarity);
   statistics.add_statistic("KKT", new_iterate.errors.KKT);
   statistics.add_statistic("FJ", new_iterate.errors.FJ);
}

bool Uno::termination_criterion(TerminationStatus current_status, int iteration) const {
   return current_status != NOT_OPTIMAL || this->max_iterations <= iteration;
}

TerminationStatus Uno::check_termination(const Problem& problem, Iterate& current_iterate, double step_norm, double objective_multiplier) const {
   TerminationStatus status = NOT_OPTIMAL;

   if (current_iterate.errors.complementarity <= this->tolerance * (double) (current_iterate.x.size() + problem.number_constraints)) {
      // feasible and KKT point
      if (current_iterate.errors.KKT <= this->tolerance * std::sqrt(current_iterate.x.size())) {
         if (current_iterate.errors.constraints <= this->tolerance * (double) current_iterate.x.size()) {
            status = KKT_POINT;
         }
      }
      // infeasible and FJ point
      else if (0 < problem.number_constraints && current_iterate.errors.FJ <= this->tolerance * std::sqrt(current_iterate.x.size())) {
         status = FJ_POINT;
      }
   }
   else if (step_norm <= this->tolerance / 100.) {
      if (current_iterate.errors.constraints <= this->tolerance * (double) current_iterate.x.size()) {
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
   std::cout << "\n";
   std::cout << "UNO v1: optimization summary\n";
   std::cout << "==============================\n";

   std::cout << "Status:\t\t\t\t";
   if (this->status == KKT_POINT) {
      std::cout << "Converged with KKT point\n";
   }
   else if (this->status == FJ_POINT) {
      std::cout << "Converged with FJ point\n";
   }
   else if (this->status == FEASIBLE_SMALL_STEP) {
      std::cout << "Converged with feasible small step\n";
   }
   else if (this->status == INFEASIBLE_SMALL_STEP) {
      std::cout << "Converged with infeasible small step\n";
   }
   else { // NOT_OPTIMAL
      std::cout << "Irregular termination\n";
   }

   std::cout << "Objective value:\t\t" << this->solution.objective << "\n";
   std::cout << "Constraint residual:\t\t" << this->solution.errors.constraints << "\n";
   std::cout << "KKT residual:\t\t\t" << this->solution.errors.KKT << "\n";
   std::cout << "FJ residual:\t\t\t" << this->solution.errors.FJ << "\n";
   std::cout << "Complementarity residual:\t" << this->solution.errors.complementarity << "\n";

   std::cout << "Feasibility measure:\t\t" << this->solution.progress.infeasibility << "\n";
   std::cout << "Optimality measure:\t\t" << this->solution.progress.objective << "\n";

   if (print_solution) {
      std::cout << "Primal solution:\t\t";
      print_vector(std::cout, this->solution.x);
      std::cout << "Lower bound multipliers:\t";
      print_vector(std::cout, this->solution.multipliers.lower_bounds);
      std::cout << "Upper bound multipliers:\t";
      print_vector(std::cout, this->solution.multipliers.upper_bounds);
      std::cout << "Constraint multipliers:\t\t";
      print_vector(std::cout, this->solution.multipliers.constraints);
   }

   std::cout << "CPU time:\t\t\t" << this->cpu_time << "s\n";
   std::cout << "Iterations:\t\t\t" << this->iteration << "\n";
   std::cout << "Objective evaluations:\t\t" << this->objective_evaluations << "\n";
   std::cout << "Constraints evaluations:\t" << this->constraint_evaluations << "\n";
   std::cout << "Jacobian evaluations:\t\t" << this->jacobian_evaluations << "\n";
   std::cout << "Hessian evaluations:\t\t" << this->hessian_evaluations << "\n";
   std::cout << "Number of subproblems solved:\t" << this->number_subproblems_solved << "\n";
}
