// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "GlobalizationMechanism.hpp"
#include "model/Model.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/Expression.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   void GlobalizationMechanism::assemble_trial_iterate(const Model& model, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double primal_step_length, double dual_step_length) {
      trial_iterate.set_number_variables(current_iterate.primals.size());
      trial_iterate.multipliers.constraints.resize(current_iterate.multipliers.constraints.size());
      trial_iterate.multipliers.lower_bounds.resize(current_iterate.multipliers.lower_bounds.size());
      trial_iterate.multipliers.upper_bounds.resize(current_iterate.multipliers.upper_bounds.size());

      // take primal step
      trial_iterate.primals = current_iterate.primals + primal_step_length * direction.primals;
      // project the trial iterate onto the bounds to avoid numerical errors
      model.project_onto_variable_bounds(trial_iterate.primals);
      // take dual step: line-search carried out only on constraint multipliers. Bound multipliers updated with full step
      trial_iterate.multipliers.constraints = current_iterate.multipliers.constraints + dual_step_length * direction.multipliers.constraints;
      trial_iterate.multipliers.lower_bounds = current_iterate.multipliers.lower_bounds + direction.multipliers.lower_bounds;
      trial_iterate.multipliers.upper_bounds = current_iterate.multipliers.upper_bounds + direction.multipliers.upper_bounds;
      trial_iterate.progress.reset();
      trial_iterate.is_objective_computed = false;
      trial_iterate.is_objective_gradient_computed = false;
      trial_iterate.are_constraints_computed = false;
      trial_iterate.is_constraint_jacobian_computed = false;
      trial_iterate.status = SolutionStatus::NOT_OPTIMAL;
   }

   void GlobalizationMechanism::set_primal_statistics(Statistics& statistics, const Model& model, const Iterate& iterate) {
      if (iterate.is_objective_computed) {
         statistics.set("objective", iterate.evaluations.objective);
      }
      if (model.is_constrained()) {
         statistics.set("primal feas", iterate.progress.infeasibility);
      }
   }

   void GlobalizationMechanism::set_dual_residuals_statistics(Statistics& statistics, const Iterate& iterate) {
      statistics.set("stationarity", iterate.residuals.stationarity);
      statistics.set("complementarity", iterate.residuals.complementarity);
   }
} // namespace