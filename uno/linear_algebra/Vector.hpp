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

Norm norm_from_string(const std::string& norm_string);

void add_vectors(const std::vector<double>& x, const std::vector<double>& y, double scaling_factor, std::vector<double>& result);

template <typename T>
void initialize_vector(std::vector<T>& x, T value) {
   for (T& xi: x) {
      xi = T(value);
   }
}

void scale(std::vector<double>& x, double scaling_factor);

template <typename T>
void copy_from(std::vector<T>& destination, const std::vector<T>& source, size_t length = std::numeric_limits<size_t>::max()) {
   length = std::min(length, std::min(source.size(), destination.size()));
   const auto source_start_position = std::cbegin(source);
   const auto source_end_position = source_start_position + length;
   const auto destination_position = std::begin(destination);
   std::copy(source_start_position, source_end_position, destination_position);
}

double norm_1(const std::vector<double>& x);
double norm_2_squared(const std::vector<double>& x);
double norm_2(const std::vector<double>& x);
double norm_inf(const std::vector<double>& x);
double norm(const std::vector<double>& x, Norm norm);

// these methods take:
// - a callback as argument. This avoids forming the vector explicitly
// - an iterable range of arbitrary type (can be Range, std::vector, etc)
template <typename RANGE>
double norm_1(const std::function<double(size_t i)>& f, const RANGE& range) {
   double norm = 0.;
   for (size_t i: range) {
      norm += std::abs(f(i));
   }
   return norm;
}

template <typename RANGE>
double norm_inf(const std::vector<double>& x, const RANGE& range) {
   double norm = 0.;
   for (size_t i: range) {
      norm = std::max(norm, std::abs(x[i]));
   }
   return norm;
}

template <typename RANGE>
double norm_inf(const std::function<double(size_t i)>& f, const RANGE& range) {
   double norm = 0.;
   for (size_t i: range) {
      norm = std::max(norm, std::abs(f(i)));
   }
   return norm;
}

template <typename RANGE>
double norm_2_squared(const std::function<double(size_t i)>& f, const RANGE& range) {
   double norm = 0.;
   for (size_t i: range) {
      double x_i = f(i);
      norm += x_i * x_i;
   }
   return norm;
}

template <typename RANGE>
double norm_2(const std::function<double(size_t i)>& f, const RANGE& range) {
   return std::sqrt(norm_2_squared(f, range));
}

template <typename RANGE>
double norm(const std::function<double(size_t i)>& f, RANGE range, Norm norm) {
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
void print_vector(std::ostream& stream, const std::vector<T>& x, size_t start = 0, size_t length = std::numeric_limits<size_t>::max(),
      const char end = '\n') {
   for (size_t i = start; i < std::min(start + length, x.size()); i++) {
      stream << x[i] << " ";
   }
   stream << end;
}

template <typename T>
void print_vector(const Level& level, const std::vector<T>& x, size_t start = 0, size_t length = std::numeric_limits<size_t>::max(),
      const char end = '\n') {
   for (size_t i = start; i < std::min(start + length, x.size()); i++) {
      level << x[i] << " ";
   }
   level << end;
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
