#ifndef UNO_VECTOR_H
#define UNO_VECTOR_H

#include <iostream>
#include <limits>
#include <map>
#include <vector>
#include <functional>
#include "tools/Logger.hpp"

enum Norm {L1_NORM = 1, L2_NORM = 2, L2_SQUARED_NORM, INF_NORM};

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

double dot(const std::vector<double>& x, const std::vector<double>& y);

double norm_1(const std::vector<double>& x);
double norm_1(const std::function<double(size_t i)>& f, size_t size);

double norm_2_squared(const std::vector<double>& x);
double norm_2_squared(const std::function<double(size_t i)>& f, size_t size);
double norm_2(const std::vector<double>& x);
double norm_2(const std::function<double(size_t i)>& f, size_t size);

double norm_inf(const std::vector<double>& x, size_t start = 0, size_t length = std::numeric_limits<size_t>::max());
double norm_inf(const std::function<double(size_t i)>& f, size_t size);

template<typename T>
double norm(const T& x, Norm norm) {
   // choose the right norm
   if (norm == INF_NORM) {
      return norm_inf(x);
   }
   else if (norm == L2_NORM) {
      return norm_2(x);
   }
   else if (norm == L2_SQUARED_NORM) {
      return norm_2_squared(x);
   }
   else if (norm == L1_NORM) {
      return norm_1(x);
   }
   throw std::out_of_range("The norm is not known");
}

double norm(const std::function<double(size_t i)>& f, size_t size, Norm norm);

template<typename T>
void print_vector(std::ostream& stream, const std::vector<T>& x, size_t start = 0, size_t length = std::numeric_limits<size_t>::max(),
      const char end = '\n') {
   for (size_t i = start; i < std::min(start + length, x.size()); i++) {
      stream << x[i] << " ";
   }
   stream << end;
}

template<typename T>
void print_vector(const Level& level, const std::vector<T>& x, size_t start = 0, size_t length = std::numeric_limits<size_t>::max(),
      const char end = '\n') {
   for (size_t i = start; i < std::min(start + length, x.size()); i++) {
      level << x[i] << " ";
   }
   level << end;
}

bool in_increasing_order(const int* array, size_t length);

#endif // UNO_VECTOR_H
