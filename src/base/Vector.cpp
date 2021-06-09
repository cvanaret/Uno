#include <cmath>
#include <cassert>
#include "Vector.hpp"

void add_vectors(const std::vector<double>& x, const std::vector<double>& y, double scaling_factor, std::vector<double>& result) {
   assert(x.size() == y.size() && "Vector.add_vectors: x and y have different sizes");
   assert(x.size() <= result.size() && "Vector.add_vectors: result is not large enough");

   for (size_t i = 0; i < x.size(); i++) {
      result[i] = x[i] + scaling_factor * y[i];
   }
}

std::vector<double> add_vectors(const std::vector<double>& x, const std::vector<double>& y, double scaling_factor) {
   std::vector<double> result(x.size());
   add_vectors(x, y, scaling_factor, result);
   return result;
}

void clear(std::vector<double>& x) {
   for (double& xi: x) {
      xi = 0.;
   }
}

void clear(SparseVector& x) {
   x.clear();
}

void scale(std::vector<double>& x, double scaling_factor) {
   for (double& xi: x) {
      xi *= scaling_factor;
   }
}
void scale(SparseVector& x, double scaling_factor) {
   for (auto& element: x) {
      element.second *= scaling_factor;
   }
}

void copy_from(std::vector<double>& destination, const std::vector<double>& source) {
   assert(destination.size() == source.size());
   for (size_t i = 0; i < source.size(); i++) {
      destination[i] = source[i];
   }
}

/* compute ||x||_1 */
double norm_1(const std::vector<double>& x) {
   double norm = 0.;
   for (double xi: x) {
      norm += std::abs(xi);
   }
   return norm;
}

double norm_1(const SparseVector& x) {
   double norm = 0.;
   for (std::pair<int, double> term: x) {
      double xi = term.second;
      norm += std::abs(xi);
   }
   return norm;
}

// https://en.wikipedia.org/wiki/Matrix_norm#Special_cases
double norm_1(const std::vector<SparseVector>& m) {
   double norm = 0.;
   for (const auto& column: m) {
      double column_norm = norm_1(column);
      norm = std::max(norm, column_norm);
   }
   return norm;
}

double norm_1(const std::function<double(int i)>& f, size_t size) {
   double norm = 0.;
   for (size_t i = 0; i < size; i++) {
      norm += std::abs(f(i));
   }
   return norm;
}

/* compute ||x||^2_2 */
double norm_2_squared(const std::vector<double>& x) {
   double norm_squared = 0.;
   for (double xi: x) {
      norm_squared += xi * xi;
   }
   return norm_squared;
}

double norm_2_squared(const SparseVector& x) {
   double norm_squared = 0.;
   for (std::pair<int, double> term: x) {
      double xi = term.second;
      norm_squared += xi * xi;
   }
   return norm_squared;
}

double norm_2_squared(const std::function<double(int i)>& f, size_t size) {
   double norm = 0.;
   for (size_t i = 0; i < size; i++) {
      double x_i = f(i);
      norm += x_i * x_i;
   }
   return norm;
}

/* compute ||x||_2 */
double norm_2(const std::vector<double>& x) {
   return std::sqrt(norm_2_squared(x));
}

double norm_2(const SparseVector& x) {
   return std::sqrt(norm_2_squared(x));
}

double norm_2(const std::function<double(int i)>& f, size_t size) {
   return std::sqrt(norm_2_squared(f, size));
}

/* compute ||x||_infty */
double norm_inf(const std::vector<double>& x, size_t length) {
   double norm = 0.;
   for (size_t i = 0; i < std::min<size_t>(length, x.size()); i++) {
      norm = std::max(norm, std::abs(x[i]));
   }
   return norm;
}

double norm_inf(const SparseVector& x) {
   double norm = 0.;
   for (std::pair<int, double> term: x) {
      double xi = term.second;
      norm = std::max(norm, std::abs(xi));
   }
   return norm;
}

// https://en.wikipedia.org/wiki/Matrix_norm#Special_cases
double norm_inf(const std::vector<SparseVector>& m) {
   // compute maximum row index
   unsigned int number_rows = 0;
   for (size_t j = 0; j < m.size(); j++) {
      number_rows = std::max(number_rows, 1 + m[j].begin()->first);
   }
   // read the matrix column-wise and fill in the row_vectors norm vector
   std::vector<double> row_vectors(number_rows);
   for (size_t j = 0; j < m.size(); j++) {
      for (const auto[i, value]: m[j]) {
         row_vectors[i] += std::abs(value);
      }
   }
   // compute the maximal component of the row_vectors vector
   double norm = 0.;
   for (double& row_vector : row_vectors) {
      norm = std::max(norm, row_vector);
   }
   return norm;
}

double norm_inf(const std::function<double(int i)>& f, size_t size) {
   double norm = 0.;
   for (size_t i = 0; i < size; i++) {
      norm = std::max(norm, std::abs(f(i)));
   }
   return norm;
}

double norm(const std::function<double(int i)>& f, size_t size, Norm norm) {
   /* choose the right norm */
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

/* dot products */

double dot(const std::vector<double>& x, const std::vector<double>& y) {
   double dot = 0.;
   for (size_t i = 0; i < std::min(x.size(), y.size()); i++) {
      dot += x[i] * y[i];
   }
   return dot;
}

double dot(const std::vector<double>& x, const SparseVector& y) {
   double dot = 0.;
   for (const auto[i, yi]: y) {
      if (i < x.size()) {
         dot += x[i] * yi;
      }
      else {
         throw std::length_error("Vector.dot: x and y have different sizes");
      }
   }
   return dot;
}

double dot(const SparseVector& x, const SparseVector& y) {
   double dot = 0.;
   for (const auto[i, xi]: x) {
      try {
         dot += xi * y.at(i);
      }
      catch (std::out_of_range&) {
      }
   }
   return dot;
}

std::string join(const std::vector<std::string>& vector, const std::string& separator) {
   std::string s;
   if (!vector.empty()) {
      s.append(vector[0]);
   }
   for (size_t i = 1; i < vector.size(); i++) {
      s.append(separator);
      s.append(vector[i]);
   }
   return s;
}
