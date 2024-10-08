// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "PrimalDualInteriorPointSystem.hpp"
#include "optimization/Iterate.hpp"
#include "reformulation/OptimizationProblem.hpp"

namespace uno {
   void PrimalDualInteriorPointSystem::right_hand_side(Vector<double>& rhs) const {
      for (size_t variable_index: Range(this->problem.number_variables)) {
         rhs[variable_index] = this->problem.variable_lower_bound(variable_index) - this->current_iterate.primals[variable_index];
      }
      // ...
   }
} // namespace