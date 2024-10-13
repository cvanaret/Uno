// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "PrimalDualInteriorPointSystem.hpp"
#include "optimization/Iterate.hpp"
#include "reformulation/OptimizationProblem.hpp"
#include "symbolic/VectorView.hpp"

namespace uno {
   PrimalDualInteriorPointSystem::PrimalDualInteriorPointSystem(const OptimizationProblem& problem, Iterate& current_iterate,
         const Multipliers& current_multipliers, double barrier_parameter, bool /*use_regularization*/, double trust_region_radius, const Options&
         options):
      LagrangeNewtonSubproblem(problem, current_iterate, current_multipliers, false, trust_region_radius, options),
      number_variables(problem.number_variables + problem.number_constraints),
      barrier_parameter(barrier_parameter) { }

   void PrimalDualInteriorPointSystem::evaluate_matrix(SymmetricMatrix<size_t, double>& matrix, const WarmstartInformation& warmstart_information) const {
      // barrier Lagrangian Hessian
      if (warmstart_information.objective_changed || warmstart_information.constraints_changed) {
         // original Lagrangian Hessian
         matrix.set_dimension(problem.number_variables + problem.number_constraints);
         this->hessian_model->evaluate(this->problem, this->current_iterate.primals, current_multipliers.constraints, matrix);

         // diagonal barrier terms (grouped by variable)
         for (size_t variable_index: Range(this->problem.number_variables)) {
            double diagonal_barrier_term = 0.;
            if (is_finite(this->problem.variable_lower_bound(variable_index))) { // lower bounded
               const double distance_to_bound = this->current_iterate.primals[variable_index] - this->problem.variable_lower_bound(variable_index);
               diagonal_barrier_term += current_multipliers.lower_bounds[variable_index] / distance_to_bound;
            }
            if (is_finite(this->problem.variable_upper_bound(variable_index))) { // upper bounded
               const double distance_to_bound = this->current_iterate.primals[variable_index] - this->problem.variable_upper_bound(variable_index);
               diagonal_barrier_term += current_multipliers.upper_bounds[variable_index] / distance_to_bound;
            }
            matrix.insert(diagonal_barrier_term, variable_index, variable_index);
         }
      }

      // constraints and Jacobian
      if (warmstart_information.constraints_changed) {
         // this->problem.evaluate_constraint_jacobian(this->current_iterate, this->constraint_jacobian);
      }
   }

   void PrimalDualInteriorPointSystem::evaluate_right_hand_side(Vector<double>& rhs, const WarmstartInformation& warmstart_information) const {
      for (size_t variable_index: Range(this->problem.number_variables)) {
         rhs[variable_index] = this->problem.variable_lower_bound(variable_index) - this->current_iterate.primals[variable_index];
      }
      // ...
      /*
      this->problem.evaluate_constraints(this->current_iterate, view(rhs, this->problem.number_variables, this->problem.number_variables +
         this->problem.number_constraints));
         */

      // barrier objective gradient
      if (warmstart_information.objective_changed) {
         // original objective gradient
         // problem.evaluate_objective_gradient(current_iterate, this->objective_gradient);

         // barrier terms
         for (size_t variable_index: Range(problem.number_variables)) {
            double barrier_term = 0.;
            if (is_finite(problem.variable_lower_bound(variable_index))) { // lower bounded
               barrier_term += -this->barrier_parameter/(this->current_iterate.primals[variable_index] - this->problem.variable_lower_bound(variable_index));
               // damping
               if (not is_finite(problem.variable_upper_bound(variable_index))) {
                  barrier_term += this->damping_factor * this->barrier_parameter;
               }
            }
            if (is_finite(problem.variable_upper_bound(variable_index))) { // upper bounded
               barrier_term += -this->barrier_parameter/(this->current_iterate.primals[variable_index] - this->problem.variable_upper_bound(variable_index));
               // damping
               if (not is_finite(problem.variable_lower_bound(variable_index))) {
                  barrier_term -= this->damping_factor * this->barrier_parameter;
               }
            }
            // this->objective_gradient.insert(variable_index, barrier_term);
         }
      }
   }
} // namespace