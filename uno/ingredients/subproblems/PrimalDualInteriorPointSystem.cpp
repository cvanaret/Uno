// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "PrimalDualInteriorPointSystem.hpp"
#include "optimization/Iterate.hpp"
#include "reformulation/OptimizationProblem.hpp"
#include "symbolic/VectorView.hpp"

namespace uno {
   PrimalDualInteriorPointSystem::PrimalDualInteriorPointSystem(const OptimizationProblem& problem, Iterate& current_iterate,
         const Multipliers& current_multipliers, double barrier_parameter, bool use_regularization, double trust_region_radius, const Options& options):
      LagrangeNewtonSubproblem(problem, current_iterate, current_multipliers, use_regularization, trust_region_radius, options),
      number_variables(problem.number_variables + problem.number_constraints),
      barrier_parameter(barrier_parameter) { }

   void PrimalDualInteriorPointSystem::evaluate_right_hand_side(Vector<double>& rhs) const {
      for (size_t variable_index: Range(this->problem.number_variables)) {
         rhs[variable_index] = this->problem.variable_lower_bound(variable_index) - this->current_iterate.primals[variable_index];
      }
      // ...
      /*
      this->problem.evaluate_constraints(this->current_iterate, view(rhs, this->problem.number_variables, this->problem.number_variables +
         this->problem.number_constraints));
         */
   }
} // namespace