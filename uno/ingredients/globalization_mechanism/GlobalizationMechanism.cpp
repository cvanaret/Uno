// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "GlobalizationMechanism.hpp"

GlobalizationMechanism::GlobalizationMechanism(ConstraintRelaxationStrategy& constraint_relaxation_strategy, const Options& options) :
      constraint_relaxation_strategy(constraint_relaxation_strategy),
      tolerance(options.get_double("tolerance")),
      unbounded_objective_threshold(options.get_double("unbounded_objective_threshold")) {
}

Iterate GlobalizationMechanism::assemble_trial_iterate(Iterate& current_iterate, const Direction& direction, double primal_dual_step_length,
      double bound_dual_step_length) {
   const auto take_dual_step = [&](Iterate& iterate) {
      // take dual step: line-search carried out only on constraint multipliers. Bound multipliers updated with full step length
      add_vectors(current_iterate.multipliers.constraints, direction.multipliers.constraints, primal_dual_step_length, iterate.multipliers.constraints);
      add_vectors(current_iterate.multipliers.lower_bounds, direction.multipliers.lower_bounds, bound_dual_step_length, iterate.multipliers.lower_bounds);
      add_vectors(current_iterate.multipliers.upper_bounds, direction.multipliers.upper_bounds, bound_dual_step_length, iterate.multipliers.upper_bounds);
      iterate.multipliers.objective = direction.objective_multiplier;
   };
   if (0. < direction.norm) {
      Iterate trial_iterate(current_iterate.primals.size(), direction.multipliers.constraints.size());
      // take primal step
      add_vectors(current_iterate.primals, direction.primals, primal_dual_step_length, trial_iterate.primals);
      // take dual step
      take_dual_step(trial_iterate);
      return trial_iterate;
   }
   else {
      // d = 0, no primal step to take. Take only dual step
      take_dual_step(current_iterate);
      current_iterate.progress = {INF<double>, {}, INF<double>};
      DEBUG << "Primal step is 0. The objective and constraints will not be re-evaluated.\n";
      return current_iterate;
   }
}

bool GlobalizationMechanism::terminate_with_small_step(const Model& model, const Direction& direction, Iterate& trial_iterate) const {
   // evaluate infeasibility
   trial_iterate.evaluate_constraints(model);
   trial_iterate.residuals.infeasibility = model.compute_constraint_violation(trial_iterate.evaluations.constraints, L1_NORM);

   // terminate with a feasible point
   if (trial_iterate.residuals.infeasibility <= this->tolerance) {
      trial_iterate.status = TerminationStatus::FEASIBLE_SMALL_STEP;
      return true;
   }
   else if (direction.objective_multiplier == 0.) { // terminate with an infeasible stationary point
      trial_iterate.status = TerminationStatus::INFEASIBLE_STATIONARY_POINT;
      return true;
   }
   else { // do not terminate, infeasible non stationary
      return false;
   }
}

TerminationStatus GlobalizationMechanism::check_termination(const Model& model, Iterate& current_iterate) const {
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

   if (current_iterate.evaluations.objective < this->unbounded_objective_threshold) {
      return TerminationStatus::UNBOUNDED;
   }
   else if (optimality_complementarity && primal_feasibility) {
      if (0. < current_iterate.multipliers.objective && optimality_stationarity) {
         // feasible regular stationary point
         return TerminationStatus::FEASIBLE_KKT_POINT;
      }
      else if (feasibility_stationarity && no_trivial_duals) {
         // feasible but CQ failure
         return TerminationStatus::FEASIBLE_FJ_POINT;
      }
   }
   else if (feasibility_complementarity && feasibility_stationarity) {
      // no primal feasibility, stationary point of constraint violation
      return TerminationStatus::INFEASIBLE_STATIONARY_POINT;
   }
   return TerminationStatus::NOT_OPTIMAL;
}

size_t GlobalizationMechanism::get_hessian_evaluation_count() const {
   return this->constraint_relaxation_strategy.get_hessian_evaluation_count();
}

size_t GlobalizationMechanism::get_number_subproblems_solved() const {
   return this->constraint_relaxation_strategy.get_number_subproblems_solved();
}