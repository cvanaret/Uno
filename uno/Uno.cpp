// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project root for details.

#include <cmath>
#include "Uno.hpp"
#include "optimization/Iterate.hpp"
#include "preprocessing/Preprocessing.hpp"
#include "preprocessing/Scaling.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"
#include "tools/Timer.hpp"

Uno::Uno(GlobalizationMechanism& globalization_mechanism, const Options& options) :
      globalization_mechanism(globalization_mechanism),
      tolerance(std::stod(options.at("tolerance"))),
      max_iterations(std::stoul(options.at("max_iterations"))),
      small_step_factor(std::stod(options.at("small_step_factor"))) {
}

Result Uno::solve(const Model& model, Iterate& current_iterate, const Options& options) {
   Timer timer{};
   timer.start();
   size_t major_iterations = 0;

   std::cout << "\nProblem " << model.name << '\n';
   std::cout << model.number_variables << " variables, " << model.number_constraints << " constraints\n";
   std::cout << "Problem type: " << Model::type_to_string[model.problem_type] << "\n\n";

   Statistics statistics = Uno::create_statistics(model, options);

   // use the current point to initialize the strategies and generate the initial iterate
   this->globalization_mechanism.initialize(statistics, current_iterate);

   TerminationStatus termination_status = NOT_OPTIMAL;
   try {
      // check for convergence
      while (!this->termination_criterion(termination_status, major_iterations)) {
         statistics.new_line();
         major_iterations++;
         DEBUG << "### Outer iteration " << major_iterations << '\n';

         // compute an acceptable iterate by solving a subproblem at the current point
         auto [new_iterate, direction_norm] = this->globalization_mechanism.compute_acceptable_iterate(statistics, current_iterate);
         assert(new_iterate.is_objective_computed && "The objective function was not evaluated at the new iterate.");

         // compute the status of the new iterate
         termination_status = this->check_termination(model, new_iterate, direction_norm);

         // TODO compute complementarity
         Uno::add_statistics(statistics, model, new_iterate, major_iterations);
         if (Logger::logger_level == INFO) statistics.print_current_line();

         current_iterate = std::move(new_iterate);
      }
   }
   catch (std::exception& exception) {
      ERROR << exception.what();
   }
   if (Logger::logger_level == INFO) statistics.print_footer();
   timer.stop();

   const size_t number_subproblems_solved = this->globalization_mechanism.get_number_subproblems_solved();
   const size_t hessian_evaluation_count = this->globalization_mechanism.get_hessian_evaluation_count();
   Result result = {termination_status, std::move(current_iterate), model.number_variables, model.number_constraints, major_iterations,
         timer.get_duration(), Iterate::number_eval_objective, Iterate::number_eval_constraints, Iterate::number_eval_jacobian, hessian_evaluation_count,
          number_subproblems_solved};
   return result;
}

Statistics Uno::create_statistics(const Model& model, const Options& options) {
   Statistics statistics(options);
   statistics.add_column("major", Statistics::int_width, 1);
   statistics.add_column("minor", Statistics::int_width, 2);
   statistics.add_column("step norm", Statistics::double_width, 31);
   statistics.add_column("f", Statistics::double_width, 100);
   if (model.is_constrained()) {
      statistics.add_column("||c||", Statistics::double_width, 101);
   }
   statistics.add_column("complementarity", Statistics::double_width, 104);
   statistics.add_column("stationarity", Statistics::double_width, 105);
   return statistics;
}

void Uno::add_statistics(Statistics& statistics, const Model& model, const Iterate& new_iterate, size_t major_iterations) {
   statistics.add_statistic(std::string("major"), major_iterations);
   statistics.add_statistic("f", new_iterate.original_evaluations.objective);
   if (model.is_constrained()) {
      statistics.add_statistic("||c||", new_iterate.constraint_violation);
   }
   statistics.add_statistic("complementarity", new_iterate.complementarity_error);
   statistics.add_statistic("stationarity", new_iterate.stationarity_error);
}

bool Uno::termination_criterion(TerminationStatus current_status, size_t iteration) const {
   return current_status != NOT_OPTIMAL || this->max_iterations <= iteration;
}

TerminationStatus Uno::check_termination(const Model& model, Iterate& current_iterate, double step_norm) const {
   const size_t number_variables = current_iterate.primals.size();

   // TODO: criterion for CQ failure
   current_iterate.evaluate_constraints(model);
   current_iterate.complementarity_error = model.compute_complementarity_error(current_iterate.primals, current_iterate.original_evaluations.constraints,
         current_iterate.multipliers.constraints, current_iterate.multipliers.lower_bounds, current_iterate.multipliers.upper_bounds);

   if (current_iterate.complementarity_error <= this->tolerance * static_cast<double>(number_variables + model.number_constraints)) {
      // feasible and KKT point
      if (current_iterate.stationarity_error <= this->tolerance * std::sqrt(number_variables) &&
          current_iterate.constraint_violation <= this->tolerance * static_cast<double>(number_variables)) {
         return KKT_POINT;
      }
      // infeasible and FJ point
      else if (model.is_constrained() && current_iterate.multipliers.objective == 0. &&
         current_iterate.stationarity_error <= this->tolerance * std::sqrt(number_variables)) {
         return FJ_POINT;
      }
   }
   if (step_norm <= this->tolerance / this->small_step_factor) {
      if (current_iterate.constraint_violation <= this->tolerance * static_cast<double>(number_variables)) {
         return FEASIBLE_SMALL_STEP;
      }
      else {
         return INFEASIBLE_SMALL_STEP;
      }
   }
   return NOT_OPTIMAL;
}

void Uno::postsolve_solution(const Model& model, const Scaling& scaling, Iterate& current_iterate, TerminationStatus termination_status) {
   // remove auxiliary variables
   current_iterate.set_number_variables(model.number_variables);

   // objective value
   current_iterate.original_evaluations.objective /= scaling.get_objective_scaling();

   // unscale the multipliers and the function values
   const bool is_feasible = (termination_status == KKT_POINT || termination_status == FEASIBLE_SMALL_STEP);
   const double scaled_objective_multiplier = scaling.get_objective_scaling()*(is_feasible ? current_iterate.multipliers.objective : 1.);
   if (scaled_objective_multiplier != 0.) {
      for (size_t j = 0; j < model.number_constraints; j++) {
         current_iterate.multipliers.constraints[j] *= scaling.get_constraint_scaling(j)/scaled_objective_multiplier;
      }
      for (size_t i = 0; i < model.number_variables; i++) {
         current_iterate.multipliers.lower_bounds[i] /= scaled_objective_multiplier;
         current_iterate.multipliers.upper_bounds[i] /= scaled_objective_multiplier;
      }
   }
}

void Result::print(bool print_solution) const {
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

   std::cout << "Objective value:\t\t" << this->solution.original_evaluations.objective << '\n';
   std::cout << "Constraint violation:\t\t" << this->solution.constraint_violation << '\n';
   std::cout << "Stationarity error:\t\t" << this->solution.stationarity_error << '\n';
   std::cout << "Complementarity error:\t\t" << this->solution.complementarity_error << '\n';

   std::cout << "Infeasibility measure:\t\t" << this->solution.nonlinear_progress.infeasibility << '\n';
   std::cout << "Optimality measure:\t\t" << this->solution.nonlinear_progress.optimality << '\n';

   if (print_solution) {
      std::cout << "Primal solution:\t\t";
      print_vector(std::cout, this->solution.primals);
      std::cout << "Lower bound multipliers:\t";
      print_vector(std::cout, this->solution.multipliers.lower_bounds);
      std::cout << "Upper bound multipliers:\t";
      print_vector(std::cout, this->solution.multipliers.upper_bounds);
      if (!this->solution.multipliers.constraints.empty()) {
         std::cout << "Constraint multipliers:\t\t";
         print_vector(std::cout, this->solution.multipliers.constraints);
      }
      std::cout << "Objective multiplier:\t\t" << this->solution.multipliers.objective << '\n';
   }

   std::cout << "CPU time:\t\t\t" << this->cpu_time << "s\n";
   std::cout << "Iterations:\t\t\t" << this->iteration << '\n';
   std::cout << "Objective evaluations:\t\t" << this->objective_evaluations << '\n';
   std::cout << "Constraints evaluations:\t" << this->constraint_evaluations << '\n';
   std::cout << "Jacobian evaluations:\t\t" << this->jacobian_evaluations << '\n';
   std::cout << "Hessian evaluations:\t\t" << this->hessian_evaluations << '\n';
   std::cout << "Number of subproblems solved:\t" << this->number_subproblems_solved << '\n';
}
