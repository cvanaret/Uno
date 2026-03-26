// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "LPSolver.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/VectorView.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/Multipliers.hpp"

namespace uno {
   // because of the way we form LPs/QPs, we get the new *multipliers* back from the solver. To get the dual displacements/direction,
   // we need to subtract the current multipliers
   void LPSolver::compute_dual_displacements(const Subproblem& subproblem, Multipliers& direction_multipliers) {
      view(direction_multipliers.constraints, 0, subproblem.number_constraints) -=
         view(subproblem.current_iterate.multipliers.constraints, 0, subproblem.number_constraints);
      view(direction_multipliers.lower_bounds, 0, subproblem.number_variables) -=
         view(subproblem.current_iterate.multipliers.lower_bounds, 0, subproblem.number_variables);
      view(direction_multipliers.upper_bounds, 0, subproblem.number_variables) -=
         view(subproblem.current_iterate.multipliers.upper_bounds, 0, subproblem.number_variables);
   }
} // namespace