// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INFINITY_H
#define UNO_INFINITY_H

#include <cmath>
#include <limits>

namespace uno {
   template <typename NumericalType>
   constexpr double INF = std::numeric_limits<NumericalType>::infinity();

   template <typename NumericalType>
   inline bool is_finite(NumericalType value) {
      return std::abs(value) < INF<NumericalType>;
   }
} // namespace

#endif // UNO_INFINITY_H