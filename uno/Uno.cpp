#include <cmath>
#include "Uno.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/Preprocessing.hpp"
#include "optimization/Scaling.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"
#include "tools/Timer.hpp"

Uno::Uno(GlobalizationMechanism& globalization_mechanism, const Options& options) :
      globalization_mechanism(globalization_mechanism),
      tolerance(std::stod(options.at("tolerance"))),
      max_iterations(std::stoul(options.at("max_iterations"))),
      small_step_factor(std::stod(options.at("small_step_factor"))) {
}

Result Uno::solve(const Problem& problem, Iterate& current_iterate, bool scale_functions, bool enforce_linear_constraints) {
   Timer timer{};
   timer.start();
   size_t major_iterations = 0;

   std::cout << "\nProblem " << problem.name << "\n";
   std::cout << problem.number_variables << " variables, " << problem.number_constraints << " constraints\n";
   std::cout << "Problem type: " << Problem::type_to_string[problem.problem_type] << "\n";

   // project x into the bounds
   problem.project_point_in_bounds(current_iterate.x);

   // initialize the function scaling
   Scaling scaling(problem.number_constraints, 100.);
   // function scaling
   if (scale_functions) {
      // evaluate the gradients
      current_iterate.evaluate_objective_gradient(problem, scaling);
      current_iterate.evaluate_constraints_jacobian(problem, scaling);
      scaling.compute(current_iterate.objective_gradient, current_iterate.constraints_jacobian);
   }
   // linear constraints feasible at initial point
   if (enforce_linear_constraints) {
      Preprocessing::enforce_linear_constraints(problem, scaling, current_iterate);
   }
   Statistics statistics = Uno::create_statistics();

   // use the current point to initialize the strategies and generate the initial iterate
   this->globalization_mechanism.initialize(statistics, problem, scaling, current_iterate);

   TerminationStatus termination_status = NOT_OPTIMAL;
   try {
      // check for convergence
      while (!this->termination_criterion(termination_status, major_iterations)) {
         statistics.new_line();
         major_iterations++;
         DEBUG << "\n### Outer iteration " << major_iterations << "\n";
         DEBUG << "Current iterate\n" << current_iterate << "\n";

         // compute an acceptable iterate by solving a subproblem at the current point
         auto [new_iterate, direction_norm] = this->globalization_mechanism.compute_acceptable_iterate(statistics, problem, scaling, current_iterate);

         Uno::add_statistics(statistics, new_iterate, major_iterations);
         if (Logger::logger_level == INFO) statistics.print_current_line();

         // compute the status of the new iterate
         termination_status = this->check_termination(problem, new_iterate, direction_norm);
         current_iterate = std::move(new_iterate);
      }
   }
   catch (std::exception& exception) {
      ERROR << exception.what();
   }
   if (Logger::logger_level == INFO) statistics.print_footer();
   Uno::postsolve_solution(problem, scaling, current_iterate, termination_status);
   timer.stop();

   const size_t number_subproblems_solved = this->globalization_mechanism.get_number_subproblems_solved();
   const size_t hessian_evaluation_count = this->globalization_mechanism.get_hessian_evaluation_count();
   Result result = {termination_status, std::move(current_iterate), scaling, problem.number_variables, problem.number_constraints, major_iterations,
         timer.get_duration(), Iterate::number_eval_objective, Iterate::number_eval_constraints, Iterate::number_eval_jacobian, hessian_evaluation_count,
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

void Uno::add_statistics(Statistics& statistics, const Iterate& new_iterate, size_t major_iterations) {
   statistics.add_statistic(std::string("major"), major_iterations);
   statistics.add_statistic("f", new_iterate.objective);
   statistics.add_statistic("||c||", new_iterate.errors.constraints);
   statistics.add_statistic("complementarity", new_iterate.errors.complementarity);
   statistics.add_statistic("KKT", new_iterate.errors.KKT);
   statistics.add_statistic("FJ", new_iterate.errors.FJ);
}

bool Uno::termination_criterion(TerminationStatus current_status, size_t iteration) const {
   return current_status != NOT_OPTIMAL || this->max_iterations <= iteration;
}

TerminationStatus Uno::check_termination(const Problem& problem, const Iterate& current_iterate, double step_norm) const {
   const size_t number_variables = current_iterate.x.size();

   if (current_iterate.errors.complementarity <= this->tolerance * static_cast<double>(number_variables + problem.number_constraints)) {
      // feasible and KKT point
      if (current_iterate.errors.KKT <= this->tolerance * std::sqrt(number_variables) &&
            current_iterate.errors.constraints <= this->tolerance * static_cast<double>(number_variables)) {
         return KKT_POINT;
      }
      // infeasible and FJ point
      else if (0 < problem.number_constraints && current_iterate.errors.FJ <= this->tolerance * std::sqrt(number_variables)) {
         return FJ_POINT;
      }
   }
   if (step_norm <= this->tolerance / this->small_step_factor) {
      if (current_iterate.errors.constraints <= this->tolerance * static_cast<double>(number_variables)) {
         return FEASIBLE_SMALL_STEP;
      }
      else {
         return INFEASIBLE_SMALL_STEP;
      }
   }
   return NOT_OPTIMAL;
}

void Uno::postsolve_solution(const Problem& problem, const Scaling& scaling, Iterate& current_iterate, TerminationStatus termination_status) {
   // remove auxiliary variables
   current_iterate.adjust_number_variables(problem.number_variables);

   // objective value
   current_iterate.objective /= scaling.get_objective_scaling();

   // unscale the multipliers and the function values
   const double scaled_objective_multiplier = scaling.get_objective_scaling()*(termination_status == KKT_POINT ?
         current_iterate.multipliers.objective : 1.);
   for (size_t j = 0; j < problem.number_constraints; j++) {
      current_iterate.multipliers.constraints[j] *= scaling.get_constraint_scaling(j)/scaled_objective_multiplier;
   }
   for (size_t i = 0; i < problem.number_variables; i++) {
      current_iterate.multipliers.lower_bounds[i] /= scaled_objective_multiplier;
      current_iterate.multipliers.upper_bounds[i] /= scaled_objective_multiplier;
   }
}

void Result::print(bool print_solution) const {
   std::cout << "\n";
   std::cout << "Uno: optimization summary\n";
   std::cout << Timer::get_current_date();
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
