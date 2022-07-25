// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "SparseVector.hpp"
#include <cmath>

// free functions
double norm_1(const SparseVector<double>& x) {
   double norm = 0.;
   x.for_each_value([&](double value) {
      norm += std::abs(value);
   });
   return norm;
}

double norm_inf(const SparseVector<double>& x) {
   double norm = 0.;
   x.for_each_value([&](double value) {
      norm = std::max(norm, std::abs(value));
   });
   return norm;
}

double dot(const std::vector<double>& x, const SparseVector<double>& y) {
   double dot_product = 0.;
   y.for_each([&](size_t i, double yi) {
      assert(i < x.size() && "Vector.dot: the sparse vector y is larger than the dense vector x");
      dot_product += x[i] * yi;
   });
   return dot_product;
}

double dot(const std::vector<double>& x, const SparseVector<double>& y, const std::function<bool (size_t i)>& predicate) {
   double dot_product = 0.;
   y.for_each([&](size_t i, double yi) {
      assert(i < x.size() && "Vector.dot: the sparse vector y is larger than the dense vector x");
      if (predicate(i)) {
         dot_product += x[i] * yi;
      }
   });
   return dot_product;
}

void scale(SparseVector<double>& x, double factor) {
   if (factor == 0.) {
      x.clear();
   }
   else {
      x.transform([=](double entry) {
         return factor*entry;
      });
   }
}
