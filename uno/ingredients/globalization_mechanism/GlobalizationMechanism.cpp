// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "GlobalizationMechanism.hpp"

GlobalizationMechanism::GlobalizationMechanism(ConstraintRelaxationStrategy& constraint_relaxation_strategy) :
      constraint_relaxation_strategy(constraint_relaxation_strategy) {
}

Iterate GlobalizationMechanism::assemble_trial_iterate(Iterate& current_iterate, const Direction& direction, double step_length) {
   const auto take_dual_step = [&](Iterate& iterate) {
      // take dual step: line-search carried out only on constraint multipliers. Bound multipliers updated with full step length
      add_vectors(current_iterate.multipliers.constraints, direction.multipliers.constraints, step_length, iterate.multipliers.constraints);
      add_vectors(current_iterate.multipliers.lower_bounds, direction.multipliers.lower_bounds, 1., iterate.multipliers.lower_bounds);
      add_vectors(current_iterate.multipliers.upper_bounds, direction.multipliers.upper_bounds, 1., iterate.multipliers.upper_bounds);
      iterate.multipliers.objective = direction.objective_multiplier;
   };
   if (0. < direction.norm) {
      Iterate trial_iterate(current_iterate.primals.size(), direction.multipliers.constraints.size());
      // take primal step
      add_vectors(current_iterate.primals, direction.primals, step_length, trial_iterate.primals);
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

size_t GlobalizationMechanism::get_hessian_evaluation_count() const {
   return this->constraint_relaxation_strategy.get_hessian_evaluation_count();
}

size_t GlobalizationMechanism::get_number_subproblems_solved() const {
   return this->constraint_relaxation_strategy.get_number_subproblems_solved();
}

void GlobalizationMechanism::print_warning(const char* message) {
   WARNING << YELLOW << message << RESET << '\n';
}