// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "Vector.hpp"

Norm norm_from_string(const std::string& norm_string) {
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

void add_vectors(const std::vector<double>& x, const std::vector<double>& y, double scaling_factor, std::vector<double>& result) {
   assert(x.size() <= y.size() && "Vector.add_vectors: x is larger than y");
   assert(x.size() <= result.size() && "Vector.add_vectors: result is not large enough");

   for (size_t i = 0; i < x.size(); i++) {
      result[i] = x[i] + scaling_factor * y[i];
   }
}

void scale(std::vector<double>& x, double scaling_factor) {
   for (double& xi: x) {
      xi *= scaling_factor;
   }
}

// compute ||x||_1
double norm_1(const std::vector<double>& x) {
   double norm = 0.;
   for (double xi: x) {
      norm += std::abs(xi);
   }
   return norm;
}

// compute ||x||^2_2
double norm_2_squared(const std::vector<double>& x) {
   double norm_squared = 0.;
   for (double xi: x) {
      norm_squared += xi * xi;
   }
   return norm_squared;
}

// compute ||x||_2
double norm_2(const std::vector<double>& x) {
   return std::sqrt(norm_2_squared(x));
}

// compute ||x||_infty
double norm_inf(const std::vector<double>& x) {
   double norm = 0.;
   for (size_t i = 0; i < x.size(); i++) {
      norm = std::max(norm, std::abs(x[i]));
   }
   return norm;
}

double norm(const std::vector<double>& x, Norm norm) {
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