// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include "Uno.hpp"
#include "ingredients/constraint_relaxation_strategy/ConstraintRelaxationStrategyFactory.hpp"
#include "ingredients/globalization_mechanism/GlobalizationMechanismFactory.hpp"
#include "ingredients/globalization_strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/subproblem/SubproblemFactory.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"
#include "tools/Timer.hpp"

Uno::Uno(GlobalizationMechanism& globalization_mechanism, const Options& options) :
      globalization_mechanism(globalization_mechanism),
      tolerance(options.get_double("tolerance")),
      max_iterations(options.get_unsigned_int("max_iterations")),
      terminate_with_small_step(options.get_bool("terminate_with_small_step")),
      small_step_threshold(options.get_double("small_step_threshold")) {
}

Result Uno::solve(Statistics& statistics, const Model& model, Iterate& current_iterate) {
   Timer timer{};
   timer.start();
   size_t major_iterations = 0;

   std::cout << "\nProblem " << model.name << '\n';
   std::cout << model.number_variables << " variables, " << model.number_constraints << " constraints\n";
   std::cout << "Problem type: " << type_to_string(model.problem_type) << "\n\n";

   // use the current point to initialize the strategies and generate the initial iterate
   this->globalization_mechanism.initialize(current_iterate);

   TerminationStatus termination_status = TerminationStatus::NOT_OPTIMAL;
   try {
      // check for termination
      while (not this->termination_criterion(termination_status, major_iterations)) {
         statistics.new_line();
         major_iterations++;
         DEBUG << "### Outer iteration " << major_iterations << '\n';

         // compute an acceptable iterate by solving a subproblem at the current point
         auto [next_iterate, step_norm] = this->globalization_mechanism.compute_next_iterate(statistics, model, current_iterate);

         // compute the status of the next iterate
         termination_status = this->check_termination(model, next_iterate, step_norm);
         Uno::add_statistics(statistics, model, next_iterate, major_iterations);
         if (Logger::logger_level == INFO) statistics.print_current_line();

         current_iterate = std::move(next_iterate);
      }
   }
   catch (std::exception& exception) {
      ERROR << exception.what();
   }
   Uno::postprocess_iterate(model, current_iterate, termination_status);

   if (Logger::logger_level == INFO) statistics.print_footer();
   timer.stop();

   const size_t number_subproblems_solved = this->globalization_mechanism.get_number_subproblems_solved();
   const size_t hessian_evaluation_count = this->globalization_mechanism.get_hessian_evaluation_count();
   Result result = {termination_status, std::move(current_iterate), model.number_variables, model.number_constraints, major_iterations,
         timer.get_duration(), Iterate::number_eval_objective, Iterate::number_eval_constraints, Iterate::number_eval_jacobian, hessian_evaluation_count,
          number_subproblems_solved};
   return result;
}

void Uno::add_statistics(Statistics& statistics, const Model& model, const Iterate& iterate, size_t major_iterations) {
   statistics.add_statistic(std::string("iters"), major_iterations);
   if (iterate.is_objective_computed) {
      statistics.add_statistic("objective", iterate.evaluations.objective);
   }
   else {
      statistics.add_statistic("objective", "-");
   }
   if (model.is_constrained()) {
      statistics.add_statistic("primal infeas.", iterate.residuals.infeasibility);
   }
}

bool Uno::termination_criterion(TerminationStatus current_status, size_t iteration) const {
   return current_status != NOT_OPTIMAL || this->max_iterations <= iteration;
}

TerminationStatus Uno::check_termination(const Model& model, Iterate& current_iterate, double step_norm) const {
   // evaluate termination conditions based on optimality conditions
   const bool optimality_stationarity = (current_iterate.residuals.optimality_stationarity/current_iterate.residuals.stationarity_scaling <=
         this->tolerance);
   const bool feasibility_stationarity = (current_iterate.residuals.feasibility_stationarity/current_iterate.residuals.stationarity_scaling <=
                                         this->tolerance);
   const bool optimality_complementarity = (current_iterate.residuals.optimality_complementarity / current_iterate.residuals.complementarity_scaling <= this->tolerance);
   const bool feasibility_complementarity = (current_iterate.residuals.feasibility_complementarity / current_iterate.residuals.complementarity_scaling
         <= this->tolerance);
   const bool primal_feasibility = (current_iterate.residuals.infeasibility <= this->tolerance);
   const bool no_trivial_duals = current_iterate.multipliers.not_all_zero(model.number_variables, this->tolerance);

   DEBUG << "Termination criteria:\n";
   DEBUG << "Stationarity (optimality): " << std::boolalpha << optimality_stationarity << '\n';
   DEBUG << "Stationarity (feasibility): " << std::boolalpha << feasibility_stationarity << '\n';
   DEBUG << "Complementarity (optimality): " << std::boolalpha << optimality_complementarity << '\n';
   DEBUG << "Complementarity (feasibility): " << std::boolalpha << feasibility_complementarity << '\n';
   DEBUG << "Primal feasibility: " << std::boolalpha << primal_feasibility << '\n';
   DEBUG << "Not all zero multipliers: " << std::boolalpha << no_trivial_duals << "\n\n";

   if (optimality_complementarity && primal_feasibility) {
      if (feasibility_stationarity && no_trivial_duals) {
         // feasible but CQ failure
         return FEASIBLE_FJ_POINT;
      }
      else if (0. < current_iterate.multipliers.objective && optimality_stationarity) {
         // feasible regular stationary point
         return FEASIBLE_KKT_POINT;
      }
   }
   else if (feasibility_complementarity && feasibility_stationarity) {
      // no primal feasibility, stationary point of constraint violation
     return INFEASIBLE_STATIONARY_POINT;
   }
   // stationarity & complementarity not achieved, but we can terminate with a small step
   if (this->terminate_with_small_step && step_norm <= this->small_step_threshold && primal_feasibility) {
      return FEASIBLE_SMALL_STEP;
   }
   return NOT_OPTIMAL;
}

void Uno::postprocess_iterate(const Model& model, Iterate& iterate, TerminationStatus termination_status) {
   // in case the objective was not yet evaluated, evaluate it
   iterate.evaluate_objective(model);
   model.postprocess_solution(iterate, termination_status);
   DEBUG << "Final iterate:\n" << iterate;
}

void join(const std::vector<std::string>& vector, char separator) {
   if (not vector.empty()) {
      std::cout << vector[0];
      for (size_t i: Range(1, vector.size())) {
         std::cout << separator << ' ' << vector[i];
      }
   }
}

void Uno::print_available_strategies() {
   std::cout << "Available strategies:\n";
   std::cout << "Constraint relaxation strategies: ";
   join(ConstraintRelaxationStrategyFactory::available_strategies(), ',');
   std::cout << '\n';
   std::cout << "Globalization mechanisms: ";
   join(GlobalizationMechanismFactory::available_strategies(), ',');
   std::cout << '\n';
   std::cout << "Globalization strategies: ";
   join(GlobalizationStrategyFactory::available_strategies(), ',');
   std::cout << '\n';
   std::cout << "Subproblems: ";
   join(SubproblemFactory::available_strategies(), ',');
   std::cout << '\n';
}