// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRIMALDUALINTERIORPOINTSYSTEM_H
#define UNO_PRIMALDUALINTERIORPOINTSYSTEM_H

#include "LagrangeNewtonSubproblem.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "reformulation/OptimizationProblem.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "tools/Options.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   // forward declaration
   template <typename ElementType>
   class Vector;

   class PrimalDualInteriorPointSystem: public LagrangeNewtonSubproblem {
   public:
      PrimalDualInteriorPointSystem(const OptimizationProblem& problem, Iterate& current_iterate, const Multipliers& current_multipliers,
            double barrier_parameter, bool use_regularization, double trust_region_radius, const Options& options);

      template <typename IndexType>
      void evaluate_functions(Statistics& statistics, const Multipliers& current_multipliers, SymmetricMatrix<IndexType, double>& hessian,
            const WarmstartInformation& warmstart_information);
      template <typename IndexType, typename ElementType>
      void evaluate_matrix(SymmetricMatrix<IndexType, ElementType>& matrix) const;
      void evaluate_right_hand_side(Vector<double>& rhs) const;

      const size_t number_variables;

   protected:
      const double barrier_parameter;
      const double damping_factor{1e-4}; // TODO
   };

   template <typename IndexType>
   inline void PrimalDualInteriorPointSystem::evaluate_functions(Statistics& statistics, const Multipliers& current_multipliers,
         SymmetricMatrix<IndexType, double>& hessian, const WarmstartInformation& warmstart_information) {
      // barrier Lagrangian Hessian
      if (warmstart_information.objective_changed || warmstart_information.constraints_changed) {

      }

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

      // constraints and Jacobian
      if (warmstart_information.constraints_changed) {
         //
         // this->problem.evaluate_constraint_jacobian(this->current_iterate, this->constraint_jacobian);
      }
   }

   template <typename IndexType, typename ElementType>
   inline void PrimalDualInteriorPointSystem::evaluate_matrix(SymmetricMatrix<IndexType, ElementType>& matrix) const {
      // original Lagrangian Hessian
      Options options;
      Statistics statistics(options);
      this->hessian_model->evaluate(statistics, this->problem, this->current_iterate.primals, current_multipliers.constraints, matrix, 0, 0);

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
} // namespace

#endif // UNO_PRIMALDUALINTERIORPOINTSYSTEM_H