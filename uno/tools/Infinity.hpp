// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INFINITY_H
#define UNO_INFINITY_H

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace uno {
   constexpr double Inf = std::numeric_limits<double>::infinity();

   inline bool is_finite(double value) {
      return std::abs(value) < Inf;
   }

   inline bool is_infinite(double value) {
      return std::abs(value) == Inf;
   }

   // user bounds whose magnitude reaches this threshold are infinite (CUTEst and AMPL convention, cf. IPOPT's
   // nlp_lower_bound_inf and nlp_upper_bound_inf)
   constexpr double infinite_bound_threshold = 1e20;

   inline double normalize_infinite_bound(double bound) {
      if (bound <= -infinite_bound_threshold) {
         return -Inf;
      }
      if (infinite_bound_threshold <= bound) {
         return Inf;
      }
      return bound;
   }

   template <typename Array>
   std::vector<double> normalize_infinite_bounds(const Array& bounds) {
      std::vector<double> normalized_bounds(bounds.size());
      for (size_t index = 0; index < bounds.size(); ++index) {
         normalized_bounds[index] = normalize_infinite_bound(bounds[index]);
      }
      return normalized_bounds;
   }
} // namespace

#endif // UNO_INFINITY_H