// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PUSHTOINTERIOR_H
#define UNO_PUSHTOINTERIOR_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include "tools/Infinity.hpp"

namespace uno {
   // push a value strictly inside its bounds (IPOPT's bound_push and bound_frac): the distance to a finite bound is at least
   // absolute_factor * max(1, |bound|) and, if both bounds are finite, at most relative_factor * (upper_bound - lower_bound).
   // A value already far enough from its bounds is unchanged, so that pushing is idempotent
   [[nodiscard]] inline double push_variable_to_interior(double value, double lower_bound, double upper_bound,
         double absolute_factor, double relative_factor) {
      assert(0. < absolute_factor && 0. < relative_factor && relative_factor <= 0.5);
      const bool has_lower_bound = is_finite(lower_bound);
      const bool has_upper_bound = is_finite(upper_bound);
      if (has_lower_bound && has_upper_bound) {
         const double range = upper_bound - lower_bound;
         assert(0. < range && "push_variable_to_interior: the bounds must satisfy lower_bound < upper_bound");
         const double lower_perturbation = std::min(absolute_factor * std::max(1., std::abs(lower_bound)), relative_factor * range);
         const double upper_perturbation = std::min(absolute_factor * std::max(1., std::abs(upper_bound)), relative_factor * range);
         value = std::max(value, lower_bound + lower_perturbation);
         value = std::min(value, upper_bound - upper_perturbation);
      }
      else if (has_lower_bound) {
         value = std::max(value, lower_bound + absolute_factor * std::max(1., std::abs(lower_bound)));
      }
      else if (has_upper_bound) {
         value = std::min(value, upper_bound - absolute_factor * std::max(1., std::abs(upper_bound)));
      }
      return value;
   }
} // namespace

#endif // UNO_PUSHTOINTERIOR_H
