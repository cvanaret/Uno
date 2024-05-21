// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_NORM_H
#define UNO_NORM_H

#include <iostream>
#include <limits>
#include <vector>
#include <functional>
#include <cmath>
#include <cassert>
#include "symbolic/VectorExpression.hpp"
#include "tools/Logger.hpp"
#include "symbolic/Range.hpp"
#include "symbolic/Collection.hpp"

enum class Norm {L1, L2, L2_SQUARED, INF};

inline Norm norm_from_string(const std::string& norm_string) {
   if (norm_string == "L1") {
      return Norm::L1;
   }
   else if (norm_string == "L2") {
      return Norm::L2;
   }
   else if (norm_string == "L2_squared") {
      return Norm::L2_SQUARED;
   }
   else if (norm_string == "INF") {
      return Norm::INF;
   }
   throw std::invalid_argument("The norm " + norm_string + " is not known");
}

// norms of any array with elements of any type

// compute l1 norm = sum |x|_i
template <typename ElementType, typename Indices, typename Callable>
ElementType norm_1(const VectorExpression<Indices, Callable>& expression) {
   ElementType norm{0};
   for (const auto [_, index]: expression) {
      norm += std::abs(expression[index]);
   }
   return norm;
}

template <typename Array, typename ElementType = typename Array::value_type>
ElementType norm_1(const Array& x) {
   ElementType norm{0};
   for (size_t index: Range(x.size())) {
      norm += std::abs(x[index]);
   }
   return norm;
}

// l1 norm of several arrays
template<typename Array, typename... Arrays>
typename Array::value_type norm_1(const Array& x, const Arrays&... other_arrays) {
   return norm_1(x) + norm_1(other_arrays...);
}

// compute l2 squared norm = sum x_i^2
template <typename Array, typename ElementType = typename Array::value_type>
ElementType norm_2_squared(const Array& x) {
   ElementType norm{0};
   for (size_t index: Range(x.size())) {
      const ElementType xi = x[index];
      norm += xi * xi;
   }
   return norm;
}

template <typename ElementType, typename Indices, typename Callable>
ElementType norm_2_squared(const VectorExpression<Indices, Callable>& expression) {
   ElementType norm{0};
   for (const auto [_, index]: expression) {
      const ElementType xi = expression[index];
      norm += xi * xi;
   }
   return norm;
}

// l2 squared norm of several arrays
template<typename Array, typename... Arrays>
typename Array::value_type norm_2_squared(const Array& x, const Arrays&... other_arrays) {
   return norm_2_squared(x) + norm_2_squared(other_arrays...);
}

// compute ||x||_2
template <typename Array>
typename Array::value_type norm_2(const Array& x) {
   return std::sqrt(norm_2_squared(x));
}

// l2 norm of several arrays
template<typename Array, typename... Arrays>
typename Array::value_type norm_2(const Array& x, const Arrays&... other_arrays) {
   return std::sqrt(norm_2_squared(x) + norm_2_squared(other_arrays...));
}

// compute ||x||_inf
template <typename ElementType, typename Indices, typename Callable>
ElementType norm_inf(const VectorExpression<Indices, Callable>& expression) {
   ElementType norm{0};
   for (const auto [_, index]: expression) {
      norm = std::max(norm, std::abs(expression[index]));
   }
   return norm;
}

template <typename Array, typename ElementType = typename Array::value_type>
ElementType norm_inf(const Array& x) {
   ElementType norm{0};
   for (size_t index: Range(x.size())) {
      norm = std::max(norm, std::abs(x[index]));
   }
   return norm;
}

template <typename Array, RangeDirection Direction, typename ElementType = typename Array::value_type>
ElementType norm_inf(const Array& x, const Range<Direction>& range) {
   ElementType norm{0};
   for (size_t index: range) {
      norm = std::max(norm, std::abs(x[index]));
   }
   return norm;
}

// inf norm of several arrays
template<typename Array, typename... Arrays>
typename Array::value_type norm_inf(const Array& x, const Arrays&... other_arrays) {
   return std::max(norm_inf(x), norm_inf(other_arrays...));
}

// norm of at least one array
template<typename Array, typename... Arrays>
typename Array::value_type norm(Norm norm, const Array& x, const Arrays&... other_arrays) {
   // choose the right norm
   if (norm == Norm::L1) {
      return norm_1(x, other_arrays...);
   }
   else if (norm == Norm::L2) {
      return norm_2(x, other_arrays...);
   }
   else if (norm == Norm::L2_SQUARED) {
      return norm_2_squared(x, other_arrays...);
   }
   else if (norm == Norm::INF) {
      return norm_inf(x, other_arrays...);
   }
   throw std::invalid_argument("The norm is not known");
}

#endif // UNO_NORM_H