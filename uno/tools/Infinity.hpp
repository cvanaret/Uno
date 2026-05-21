// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INFINITY_H
#define UNO_INFINITY_H

#include <cmath>
#include <limits>

namespace uno {
   template <typename T>
   const double INF = std::numeric_limits<T>::infinity();

   template <typename T>
   bool is_finite(T value) {
      return std::abs(value) < INF<T>;
   }

   template <typename T>
   bool is_infinite(T value) {
      return std::abs(value) == INF<T>;
   }
} // namespace

#endif // UNO_INFINITY_H