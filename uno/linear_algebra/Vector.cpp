#include <cmath>
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

double dot(const std::vector<double>& x, const std::vector<double>& y) {
   assert(x.size() == y.size() && "Vector.dot: x and y have different sizes");
   double result = 0.;
   for (size_t i = 0; i < x.size(); i++) {
      result += x[i]*y[i];
   }
   return result;
}

// compute ||x||_1
double norm_1(const std::vector<double>& x) {
   double norm = 0.;
   for (double xi: x) {
      norm += std::abs(xi);
   }
   return norm;
}

// this version takes a callback as argument. This avoids forming the vector explicitly
double norm_1(const std::function<double(size_t i)>& f, size_t size) {
   double norm = 0.;
   for (size_t i = 0; i < size; i++) {
      norm += std::abs(f(i));
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

// this version takes a callback as argument. This avoids forming the vector explicitly
double norm_2_squared(const std::function<double(size_t i)>& f, size_t size) {
   double norm = 0.;
   for (size_t i = 0; i < size; i++) {
      double x_i = f(i);
      norm += x_i * x_i;
   }
   return norm;
}

// compute ||x||_2
double norm_2(const std::vector<double>& x) {
   return std::sqrt(norm_2_squared(x));
}

double norm_2(const std::function<double(size_t i)>& f, size_t size) {
   return std::sqrt(norm_2_squared(f, size));
}

// compute ||x||_infty
double norm_inf(const std::vector<double>& x, size_t start, size_t length) {
   double norm = 0.;
   for (size_t i = start; i < std::min<size_t>(start + length, x.size()); i++) {
      norm = std::max(norm, std::abs(x[i]));
   }
   return norm;
}

double norm_inf(const std::function<double(size_t i)>& f, size_t size) {
   double norm = 0.;
   for (size_t i = 0; i < size; i++) {
      norm = std::max(norm, std::abs(f(i)));
   }
   return norm;
}

double norm(const std::function<double(size_t i)>& f, size_t size, Norm norm) {
   // choose the right norm
   if (norm == INF_NORM) {
      return norm_inf(f, size);
   }
   else if (norm == L2_NORM) {
      return norm_2(f, size);
   }
   else if (norm == L2_SQUARED_NORM) {
      return norm_2_squared(f, size);
   }
   else if (norm == L1_NORM) {
      return norm_1(f, size);
   }
   else {
      throw std::out_of_range("The norm is not known");
   }
}

// check that an array of integers is in increasing order (x[i] <= x[i+1])
bool in_increasing_order(const int* array, size_t length) {
   size_t i = 0;
   while (i < length-1) {
      if (array[i] > array[i+1]) {
         return false;
      }
      i++;
   }
   return true;
}