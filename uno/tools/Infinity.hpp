// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INFINITY_H
#define UNO_INFINITY_H

#include <cmath>
#include <limits>

namespace uno {
   constexpr double Inf = std::numeric_limits<double>::infinity();

   inline bool is_finite(double value) {
      return std::abs(value) < Inf;
   }

   inline bool is_infinite(double value) {
      return std::abs(value) == Inf;
   }
} // namespace

#endif // UNO_INFINITY_H