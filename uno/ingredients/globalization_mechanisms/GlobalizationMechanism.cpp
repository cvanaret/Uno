// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "GlobalizationMechanism.hpp"
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategy.hpp"
#include "model/Model.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/Expression.hpp"
#include "options/Options.hpp"

namespace uno {
   GlobalizationMechanism::GlobalizationMechanism(ConstraintRelaxationStrategy& constraint_relaxation_strategy) :
         constraint_relaxation_strategy(constraint_relaxation_strategy),
         direction(this->constraint_relaxation_strategy.maximum_number_variables(), this->constraint_relaxation_strategy.maximum_number_constraints()) {
   }

   void GlobalizationMechanism::assemble_trial_iterate(const Model& model, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double primal_step_length, double dual_step_length) {
      trial_iterate.set_number_variables(current_iterate.primals.size());
      // take primal step
      trial_iterate.primals = current_iterate.primals + primal_step_length * direction.primals;
      // project the trial iterate onto the bounds to avoid numerical errors
      model.project_onto_variable_bounds(trial_iterate.primals);
      // take dual step: line-search carried out only on constraint multipliers. Bound multipliers updated with full step
      trial_iterate.multipliers.constraints = current_iterate.multipliers.constraints + dual_step_length * direction.multipliers.constraints;
      trial_iterate.multipliers.lower_bounds = current_iterate.multipliers.lower_bounds + direction.multipliers.lower_bounds;
      trial_iterate.multipliers.upper_bounds = current_iterate.multipliers.upper_bounds + direction.multipliers.upper_bounds;
      trial_iterate.feasibility_multipliers.constraints = current_iterate.feasibility_multipliers.constraints + dual_step_length * direction.feasibility_multipliers.constraints;
      trial_iterate.feasibility_multipliers.lower_bounds = current_iterate.feasibility_multipliers.lower_bounds + direction.feasibility_multipliers.lower_bounds;
      trial_iterate.feasibility_multipliers.upper_bounds = current_iterate.feasibility_multipliers.upper_bounds + direction.feasibility_multipliers.upper_bounds;
      trial_iterate.progress.reset();
      trial_iterate.is_objective_computed = false;
      trial_iterate.is_objective_gradient_computed = false;
      trial_iterate.are_constraints_computed = false;
      trial_iterate.is_constraint_jacobian_computed = false;
      trial_iterate.status = IterateStatus::NOT_OPTIMAL;
   }

   size_t GlobalizationMechanism::get_hessian_evaluation_count() const {
      return this->constraint_relaxation_strategy.get_hessian_evaluation_count();
   }

   size_t GlobalizationMechanism::get_number_subproblems_solved() const {
      return this->constraint_relaxation_strategy.get_number_subproblems_solved();
   }
} // namespace
