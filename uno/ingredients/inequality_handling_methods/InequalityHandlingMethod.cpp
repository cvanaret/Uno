#include "InequalityHandlingMethod.hpp"
#include "optimization/Multipliers.hpp"
#include "symbolic/VectorView.hpp"

namespace uno {
   void InequalityHandlingMethod::compute_dual_displacements(const Multipliers& current_multipliers, Multipliers& direction_multipliers) {
      // compute dual *displacements* (active-set methods usually compute the new duals, not the displacements)
      view(direction_multipliers.constraints, 0, current_multipliers.constraints.size()) -= current_multipliers.constraints;
      view(direction_multipliers.lower_bounds, 0, current_multipliers.lower_bounds.size()) -= current_multipliers.lower_bounds;
      view(direction_multipliers.upper_bounds, 0, current_multipliers.upper_bounds.size()) -= current_multipliers.upper_bounds;
   }
} // namespace