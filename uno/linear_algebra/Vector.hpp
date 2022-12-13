// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VECTOR_H
#define UNO_VECTOR_H

#include <iostream>
#include <limits>
#include <vector>
#include <functional>
#include <cmath>
#include "tools/Logger.hpp"

enum Norm {
   L1_NORM = 1, L2_NORM = 2, L2_SQUARED_NORM, INF_NORM
};

inline Norm norm_from_string(const std::string& norm_string) {
   if (norm_string == "L1") {
      return L1_NORM;
   }
   else if (norm_string == "L2") {
      return L2_NORM;
   }
   else if (norm_string == "INF") {
      return INF_NORM;
   }
   throw std::out_of_range("The norm is not known");
}

// result <- x + scaling_factor * y
template <typename T>
void add_vectors(const std::vector<T>& x, const std::vector<T>& y, T scaling_factor, std::vector<T>& result) {
   assert(x.size() <= y.size() && "Vector.add_vectors: x is longer than y");
   assert(x.size() <= result.size() && "Vector.add_vectors: result is not long enough");

   for (size_t i = 0; i < x.size(); i++) {
      result[i] = x[i] + scaling_factor * y[i];
   }
}

template <typename T>
void initialize_vector(std::vector<T>& x, T value) {
   for (T& xi: x) {
      xi = T(value);
   }
}

template <typename T>
void scale(std::vector<T>& x, T scaling_factor) {
   for (T& xi: x) {
      xi *= scaling_factor;
   }
}

template <typename T>
T dot(const std::vector<T>& x, const std::vector<T>& y) {
   assert(x.size() == y.size() && "The vectors do not have the same size.");

   T dot_product = 0.;
   for (size_t i = 0; i < x.size(); i++) {
      dot_product += x[i]*y[i];
   }
   return dot_product;
}

template <typename T>
void copy_from(std::vector<T>& destination, const std::vector<T>& source, size_t length = std::numeric_limits<size_t>::max()) {
   length = std::min(length, std::min(source.size(), destination.size()));
   const auto source_start_position = std::cbegin(source);
   const auto source_end_position = source_start_position + length;
   const auto destination_position = std::begin(destination);
   std::copy(source_start_position, source_end_position, destination_position);
}

// compute ||x||_1
template <typename T>
T norm_1(const std::vector<T>& x) {
   T norm = T(0);
   for (T xi: x) {
      norm += std::abs(xi);
   }
   return norm;
}

// compute ||x||^2_2
template <typename T>
T norm_2_squared(const std::vector<T>& x) {
   T norm_squared = T(0);
   for (T xi: x) {
      norm_squared += xi * xi;
   }
   return norm_squared;
}

// compute ||x||_2
template <typename T>
T norm_2(const std::vector<T>& x) {
   return std::sqrt(norm_2_squared(x));
}

// compute ||x||_inf
template <typename T>
T norm_inf(const std::vector<T>& x) {
   T norm = T(0);
   for (T xi: x) {
      norm = std::max(norm, std::abs(xi));
   }
   return norm;
}

template <typename T>
T norm(const std::vector<T>& x, Norm norm) {
   // choose the right norm
   if (norm == INF_NORM) {
      return norm_inf<T>(x);
   }
   else if (norm == L2_NORM) {
      return norm_2<T>(x);
   }
   else if (norm == L2_SQUARED_NORM) {
      return norm_2_squared<T>(x);
   }
   else if (norm == L1_NORM) {
      return norm_1<T>(x);
   }
   throw std::out_of_range("The norm is not known");
}

// these methods take:
// - a callback as argument whose parameter is the current index. This avoids forming the vector explicitly
// - an iterable range of arbitrary type (can be Range, std::vector, etc)
template <typename T, typename RANGE>
T norm_1(const std::function<T(size_t i)>& f, const RANGE& range, bool normalized = false) {
   T norm = T(0);
   for (size_t i: range) {
      norm += std::abs(f(i));
   }
   // if required, normalize with the size of the range
   if (normalized && 0 < range.size()) {
      norm /= static_cast<double>(range.size());
   }
   return norm;
}

template <typename T, typename RANGE>
T norm_inf(const std::vector<T>& x, const RANGE& range, bool /*normalized*/ = false) {
   T norm = T(0);
   for (size_t i: range) {
      norm = std::max(norm, std::abs(x[i]));
   }
   // already normalized
   return norm;
}

template <typename T, typename RANGE>
T norm_inf(const std::function<T(size_t i)>& f, const RANGE& range, bool /*normalized*/ = false) {
   T norm = T(0);
   for (size_t i: range) {
      norm = std::max(norm, std::abs(f(i)));
   }
   // already normalized
   return norm;
}

template <typename T, typename RANGE>
T norm_2_squared(const std::function<T(size_t i)>& f, const RANGE& range, bool normalized = false) {
   T norm = T(0);
   for (size_t i: range) {
      const T x_i = f(i);
      norm += x_i * x_i;
   }
   // if required, normalize with the size of the range
   if (normalized && 0 < range.size()) {
      norm /= static_cast<double>(range.size());
   }
   return norm;
}

template <typename T, typename RANGE>
T norm_2(const std::function<T(size_t /*i*/)>& f, const RANGE& range, bool normalized = false) {
   return std::sqrt(norm_2_squared(f, range, normalized));
}

template <typename T, typename RANGE>
T norm(const std::function<T(size_t /*i*/)>& f, RANGE range, Norm norm) {
   // choose the right norm
   if (norm == INF_NORM) {
      return norm_inf(f, range);
   }
   else if (norm == L2_NORM) {
      return norm_2(f, range);
   }
   else if (norm == L2_SQUARED_NORM) {
      return norm_2_squared(f, range);
   }
   else if (norm == L1_NORM) {
      return norm_1(f, range);
   }
   else {
      throw std::out_of_range("The norm is not known");
   }
}

template <typename T>
void print_vector(std::ostream& stream, const std::vector<T>& x, size_t start = 0, size_t length = std::numeric_limits<size_t>::max()) {
   for (size_t i = start; i < std::min(start + length, x.size()); i++) {
      stream << x[i] << " ";
   }
   stream << '\n';
}

template <typename T>
void print_vector(const Level& level, const std::vector<T>& x, size_t start = 0, size_t length = std::numeric_limits<size_t>::max()) {
   for (size_t i = start; i < std::min(start + length, x.size()); i++) {
      level << x[i] << " ";
   }
   level << '\n';
}

// check that an array of integers is in increasing order (x[i] <= x[i+1])
template <typename ARRAY>
bool in_increasing_order(const ARRAY& array, size_t length) {
   size_t i = 0;
   while (i < length-1) {
      if (array[i] > array[i+1]) {
         return false;
      }
      i++;
   }
   return true;
}

#endif // UNO_VECTOR_H
