// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VECTOR_H
#define UNO_VECTOR_H

#include <iostream>
#include <limits>
#include <vector>
#include <functional>
#include <cmath>
#include "tools/Logger.hpp"
#include "tools/Range.hpp"

enum class Norm {L1, L2, L2_SQUARED, INF};

inline Norm norm_from_string(const std::string& norm_string) {
   if (norm_string == "L1") {
      return Norm::L1;
   }
   else if (norm_string == "L2") {
      return Norm::L2;
   }
   else if (norm_string == "INF") {
      return Norm::INF;
   }
   throw std::invalid_argument("The norm " + norm_string + " is not known");
}

// result <- x + scaling_factor * y
template <typename ElementType>
void add_vectors(const std::vector<ElementType>& x, const std::vector<ElementType>& y, ElementType scaling_factor, std::vector<ElementType>& result) {
   assert(x.size() <= y.size() && "Vector.add_vectors: x is longer than y");
   assert(x.size() <= result.size() && "Vector.add_vectors: result is not long enough");

   for (size_t index: Range(x.size())) {
      result[index] = x[index] + scaling_factor * y[index];
   }
}

// result <- a*x + b*y
template <typename ElementType>
void add_vectors(ElementType a, const std::vector<ElementType>& x, ElementType b, const std::vector<ElementType>& y, std::vector<ElementType>& result) {
   assert(x.size() <= y.size() && "Vector.add_vectors: x is longer than y");
   assert(x.size() <= result.size() && "Vector.add_vectors: result is not long enough");

   for (size_t index: Range(x.size())) {
      result[index] = a*x[index] + b*y[index];
   }
}

template <typename ElementType>
void initialize_vector(std::vector<ElementType>& x, ElementType value) {
   for (ElementType& xi: x) {
      xi = ElementType(value);
   }
}

template <typename ElementType>
void scale(std::vector<ElementType>& x, ElementType scaling_factor) {
   for (ElementType& xi: x) {
      xi *= scaling_factor;
   }
}

template <typename ElementType>
ElementType dot_product(const std::vector<ElementType>& x, const std::vector<ElementType>& y) {
   assert(x.size() == y.size() && "The vectors do not have the same size.");

   ElementType result = ElementType(0);
   for (size_t index: Range(x.size())) {
      result += x[index]*y[index];
   }
   return result;
}

template <typename ElementType>
void copy_from(std::vector<ElementType>& destination, const std::vector<ElementType>& source, size_t length = std::numeric_limits<size_t>::max()) {
   length = std::min(length, std::min(source.size(), destination.size()));
   const auto source_start_position = std::cbegin(source);
   const auto source_end_position = source_start_position + length;
   const auto destination_position = std::begin(destination);
   std::copy(source_start_position, source_end_position, destination_position);
}

// norms of any array with elements of any type

// compute l1 norm = sum |x|_i
template <typename Array, typename ElementType = typename Array::value_type>
ElementType norm_1(const Array& x) {
   ElementType norm{0};
   for (size_t index = 0; index < x.size(); index++) {
      norm += std::abs(x[index]);
   }
   return norm;
}

// l1 norm of several arrays
template<typename Array, typename... Arrays, typename ElementType = typename Array::value_type>
ElementType norm_1(const Array& x, Arrays... other_arrays) {
   return norm_1(x) + norm_1(other_arrays...);
}

// compute l2 squared norm = sum x_i^2
template <typename Array, typename ElementType = typename Array::value_type>
ElementType norm_2_squared(const Array& x) {
   ElementType norm_squared{0};
   for (size_t index = 0; index < x.size(); index++) {
      const ElementType xi = x[index];
      norm_squared += xi * xi;
   }
   return norm_squared;
}

// l2 squared norm of several arrays
template<typename Array, typename... Arrays, typename ElementType = typename Array::value_type>
ElementType norm_2_squared(const Array& x, Arrays... other_arrays) {
   return norm_2_squared(x) + norm_2_squared(other_arrays...);
}

// compute ||x||_2
template <typename Array, typename ElementType = typename Array::value_type>
ElementType norm_2(const Array& x) {
   return std::sqrt(norm_2_squared(x));
}

// l2 norm of several arrays
template<typename Array, typename... Arrays, typename ElementType = typename Array::value_type>
ElementType norm_2(const Array& x, Arrays... other_arrays) {
   return std::sqrt(norm_2_squared(x) + norm_2_squared(other_arrays...));
}

// compute ||x||_inf
template <typename Array, typename ElementType = typename Array::value_type>
ElementType norm_inf(const Array& x) {
   ElementType norm{0};
   for (size_t index = 0; index < x.size(); index++) {
      norm = std::max(norm, std::abs(x[index]));
   }
   return norm;
}

// inf norm of several arrays
template<typename Array, typename... Arrays, typename ElementType = typename Array::value_type>
ElementType norm_inf(const Array& x, Arrays... other_arrays) {
   return std::max(norm_inf(x), norm_inf(other_arrays...));
}

// inf norm where the indices live in a given set
template <typename ElementType, typename Array>
ElementType norm_inf(const std::vector<ElementType>& x, const Array& indices) {
   ElementType norm = ElementType(0);
   for (size_t index: indices) {
      norm = std::max(norm, std::abs(x[index]));
   }
   return norm;
}

// norm of at least one array
template<typename Array, typename... Arrays, typename ElementType = typename Array::value_type>
ElementType norm(Norm norm, const Array& x, Arrays... other_arrays) {
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

// use && to allow temporaries (such as std::cout or logger DEBUG, WARNING, etc)
template <typename Array, typename Stream>
void print_vector(Stream&& stream, const Array& x, size_t start = 0, size_t length = std::numeric_limits<size_t>::max()) {
   for (size_t index: Range(start, std::min(start + length, x.size()))) {
      stream << x[index] << " ";
   }
   stream << '\n';
}

// check that an array of integers is in increasing order (x[i] <= x[i+1])
template <typename Array>
bool in_increasing_order(const Array& array, size_t length) {
   size_t index = 0;
   while (index < length - 1) {
      if (array[index] > array[index + 1]) {
         return false;
      }
      index++;
   }
   return true;
}

#endif // UNO_VECTOR_H
