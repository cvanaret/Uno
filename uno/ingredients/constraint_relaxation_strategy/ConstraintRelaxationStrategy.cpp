// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project root for details.

#include "ConstraintRelaxationStrategy.hpp"

ConstraintRelaxationStrategy::ConstraintRelaxationStrategy(bool penalty_parameter_control, const Options& options):
      penalty_parameter_control(penalty_parameter_control), residual_norm(norm_from_string(options.at("residual_norm"))),
      small_step_threshold(stod(options.at("small_step_threshold"))) {
}

bool ConstraintRelaxationStrategy::is_small_step(const Direction& direction) const {
   return (direction.norm <= this->small_step_threshold);
}

void ConstraintRelaxationStrategy::compute_nonlinear_residuals(const ReformulatedProblem& problem, Iterate& iterate) const {
   iterate.evaluate_constraints(problem.model);
   iterate.constraint_violation = problem.model.compute_constraint_violation(iterate.original_evaluations.constraints, L1_NORM);
   iterate.evaluate_lagrangian_gradient(problem.model, problem.get_objective_multiplier(), iterate.multipliers.constraints,
         iterate.multipliers.lower_bounds, iterate.multipliers.upper_bounds);
   iterate.stationarity_error = norm(iterate.lagrangian_gradient, this->residual_norm);
}

void ConstraintRelaxationStrategy::recover_active_set(const Model& model, Direction& direction) {
   // TODO remove elastic variables
   for (size_t i = model.number_variables; i < direction.primals.size(); i++) {
      //direction.active_set.bounds.at_lower_bound.erase(i);
      //direction.active_set.bounds.at_upper_bound.erase(i);
   }
   // constraints: only when p-n = 0
   if (direction.constraint_partition.has_value()) {
      /* TODO with this->relaxed_problem
      ConstraintPartition& constraint_partition = direction.constraint_partition.value();
      this->elastic_variables.positive.for_each([&](size_t j, size_t i) {
         // if the component is strictly positive, the constraint is violated
         if (0. < direction.x[i]) {
            constraint_partition.lower_bound_infeasible.push_back(j);
         }
      });
      this->elastic_variables.negative.for_each([&](size_t j, size_t i) {
         // if the component is strictly positive, the constraint is violated
         if (0. < direction.x[i]) {
            constraint_partition.upper_bound_infeasible.push_back(j);
         }
      });
       */
   }

   /*
   for (size_t j = 0; j < direction.multipliers.constraints.size(); j++) {
      // compute constraint violation
      double constraint_violation = 0.;
      try {
         size_t i = this->elastic_variables.positive.at(j);
         constraint_violation += direction.x[i];
      }
      catch (const std::out_of_range& e) {
      }
      try {
         size_t i = this->elastic_variables.negative.at(j);
         constraint_violation += direction.x[i];
      }
      catch (const std::out_of_range& e) {
      }
      // update active set
      if (0. < constraint_violation) {
         //direction.active_set.constraints.at_lower_bound.erase(j);
         //direction.active_set.constraints.at_upper_bound.erase(j);
      }
   }
    */
}