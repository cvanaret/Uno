// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INFINITY_H
#define UNO_INFINITY_H

#include <cmath>
#include <limits>

template <typename T>
const double INF = std::numeric_limits<T>::infinity();

template <typename T>
inline bool is_finite(T value) {
   return std::abs(value) < INF<T>;
}

#endif // UNO_INFINITY_H