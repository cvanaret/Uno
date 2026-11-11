// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INITIALPOINTPERTURBATION_H
#define UNO_INITIALPOINTPERTURBATION_H

#include <algorithm>
#include <cmath>
#include <vector>
#include "linear_algebra/Vector.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   // deterministic perturbation of a point whose functions or derivatives cannot be evaluated (e.g. a square root or a
   // distance evaluated at 0, two coincident coordinates). Each non-fixed variable is moved by
   // relative_magnitude * max(1, |x_i|), scaled by a factor in [1, 2) and with a sign that both depend on the index i (so
   // that coordinates with equal values are separated). If the move leaves the bounds, the opposite direction is taken;
   // the result is projected onto the bounds
   inline void perturb_point(Vector<double>& x, const std::vector<double>& lower_bounds, const std::vector<double>& upper_bounds,
         double relative_magnitude) {
      // golden-ratio sequence: equidistributed in [0, 1), deterministic and free of short cycles
      constexpr double golden_ratio_conjugate = 0.6180339887498949;
      for (size_t variable_index: Range(x.size())) {
         const double lower_bound = lower_bounds[variable_index];
         const double upper_bound = upper_bounds[variable_index];
         if (lower_bound == upper_bound) {
            continue; // fixed variable
         }
         const double fraction = std::fmod(static_cast<double>(variable_index + 1) * golden_ratio_conjugate, 1.);
         const double sign = (fraction < 0.5) ? 1. : -1.;
         const double distance = relative_magnitude * std::max(1., std::abs(x[variable_index])) * (1. + fraction);
         double perturbed_value = x[variable_index] + sign * distance;
         if (perturbed_value < lower_bound || upper_bound < perturbed_value) {
            perturbed_value = x[variable_index] - sign * distance;
         }
         x[variable_index] = std::min(std::max(perturbed_value, lower_bound), upper_bound);
      }
   }
} // namespace

#endif // UNO_INITIALPOINTPERTURBATION_H
