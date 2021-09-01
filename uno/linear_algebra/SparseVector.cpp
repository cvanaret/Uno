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

double dot(const std::vector<double>& x, const SparseVector<double>& y) {
   double dot_product = 0.;
   y.for_each([&](size_t i, double yi) {
      assert(i < x.size() && "Vector.dot: x and y have different sizes");
      dot_product += x[i] * yi;
   });
   return dot_product;
}

void scale(SparseVector<double>& x, double factor) {
   x.transform([=](double entry) {
      return factor*entry;
   });
}