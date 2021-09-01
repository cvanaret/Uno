#include "SparseVector.hpp"
#include <cmath>

// free functions

double norm_1(const SparseVector2<double>& x) {
   double norm = 0.;
   x.for_each_value([&](double value) {
      norm += std::abs(value);
   });
   return norm;
}

double dot(const std::vector<double>& x, const SparseVector2<double>& y) {
   double dot = 0.;
   y.for_each([&](size_t i, double yi) {
      if (i < x.size()) {
         dot += x[i] * yi;
      }
      else {
         throw std::length_error("Vector.dot: x and y have different sizes");
      }
   });
   return dot;
}

void scale(SparseVector2<double>& x, double factor) {
   x.transform([=](double entry) {
      return factor*entry;
   });
}