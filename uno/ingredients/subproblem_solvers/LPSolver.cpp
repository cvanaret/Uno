// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "LPSolver.hpp"
#include "linear_algebra/VectorView.hpp"
#include "optimization/Multipliers.hpp"

namespace uno {
   void LPSolver::compute_dual_displacements(const Multipliers& current_multipliers, Multipliers& direction_multipliers) {
      view(direction_multipliers.constraints, 0, current_multipliers.constraints.size()) -= current_multipliers.constraints;
      view(direction_multipliers.lower_bounds, 0, current_multipliers.lower_bounds.size()) -= current_multipliers.lower_bounds;
      view(direction_multipliers.upper_bounds, 0, current_multipliers.upper_bounds.size()) -= current_multipliers.upper_bounds;
   }
} // namespace