#include <exception>
#include <cmath>
#include <cassert>
#include <iomanip>
#include "SymmetricMatrix.hpp"
#include "Vector.hpp"

SymmetricMatrix::SymmetricMatrix(size_t dimension, size_t capacity) : dimension(dimension), capacity(capacity) {
}

void SymmetricMatrix::reset() {
   this->number_nonzeros = 0;
}

std::vector<double> SymmetricMatrix::product(const std::vector<double>& vector) const {
   assert(this->dimension == (int) vector.size() && "The matrix and the vector do not have the same size");

   /* create (n, 1) result */
   std::vector<double> result(vector.size());
   this->for_each([&](size_t i, size_t j, double entry) {
      result[i] += entry * vector[j];
      // off-diagonal term
      if (i != j) {
         result[j] += entry * vector[i];
      }
   });
   return result;
}

double SymmetricMatrix::quadratic_product(const std::vector<double>& x, const std::vector<double>& y) const {
   assert(x.size() == y.size() && "The two vectors x and y do not have the same size");

   double result = 0.;
   this->for_each([&](size_t i, size_t j, double entry) {
      result += (i == j ? 1 : 2) * entry * x[i] * y[j];
   });
   return result;
}

void SymmetricMatrix::add_outer_product(const SparseVector& x, double scaling_factor) {
   /* perform matrix addition: A + rho x x^T */
   for (const auto[row_index, row_term]: x) {
      for (const auto[column_index, column_term]: x) {
         // upper triangular matrix
         if (row_index <= column_index) {
            // add product of components
            this->insert(scaling_factor * row_term * column_term, (int) row_index, (int) column_index);
         }
      }
   }
}

double SymmetricMatrix::smallest_diagonal_entry() const {
   double smallest_entry = INFINITY;
   this->for_each([&](size_t i, size_t j, double entry) {
      if (i == j) {
         smallest_entry = std::min(smallest_entry, entry);
      }
   });
   if (smallest_entry == INFINITY) {
      smallest_entry = 0.;
   }
   return smallest_entry;
}

/* 
 * Uno matrix: bijection between indices and single key + sparse vector (map)
 */

//UnoMatrix::UnoMatrix(size_t dimension, size_t capacity): SymmetricMatrix(dimension, capacity) {
//}
//
//void UnoMatrix::insert(double term, size_t row_index, size_t column_index) {
//   // generate the unique key
//   size_t key = column_index * this->dimension + row_index;
//   // insert the element
//   this->matrix[key] += term;
//   this->number_nonzeros++;
//}
//
//double UnoMatrix::norm_1() {
//   // compute maximum column index
//   size_t number_columns = 0;
//   for (std::pair<const size_t, double> element: this->matrix) {
//      size_t key = element.first;
//      // retrieve indices
//      size_t j = key / this->dimension;
//      number_columns = std::max(number_columns, 1 + j);
//   }
//   // read the matrix and fill in the column_vectors norm vector
//   std::vector<double> column_vectors(number_columns);
//   for (const auto[key, value]: this->matrix) {
//      // retrieve indices
//      size_t j = key / this->dimension;
//      column_vectors[j] += std::abs(value);
//   }
//   // compute the maximal component of the column_vectors vector
//   double norm = 0.;
//   for (double column_vector : column_vectors) {
//      norm = std::max(norm, column_vector);
//   }
//   return norm;
//}
//
//std::vector<double> UnoMatrix::product(const std::vector<double>& /*vector*/) {
//   throw std::out_of_range("UnoMatrix::product is not implemented");
//}
//
//void UnoMatrix::add_matrix(UnoMatrix& other_matrix, double factor) {
//   for (const auto[key, value]: other_matrix.matrix) {
//      // retrieve indices
//      size_t i = key % this->dimension;
//      size_t j = key / this->dimension;
//      this->insert(factor * value, i, j);
//   }
//}
//
//COOSymmetricMatrix UnoMatrix::to_COO() {
//   COOSymmetricMatrix coo_matrix(this->dimension, this->capacity);
//
//   for (const auto[key, value]: this->matrix) {
//      // retrieve indices
//      size_t i = key % this->dimension;
//      size_t j = key / this->dimension;
//      coo_matrix.insert(value, i, j);
//   }
//   return coo_matrix;
//}
//
///* generate a COO matrix by removing some variables (e.g. reduced Hessian in EQP problems) */
///* mask contains (i_origin, i_reduced) pairs, where i_origin is the original index, and i_reduced is the index in the reduced matrix */
//COOSymmetricMatrix UnoMatrix::to_COO(const std::unordered_map<size_t, size_t>& mask) {
//   COOSymmetricMatrix coo_matrix(this->dimension, this->capacity);
//
//   for (const auto[key, value]: this->matrix) {
//      // retrieve indices
//      size_t i = key % this->dimension;
//      size_t j = key / this->dimension;
//      try {
//         // if i and j are kept, compute the indices in the reduced matrix
//         size_t i_reduced = mask.at(i);
//         size_t j_reduced = mask.at(j);
//         coo_matrix.insert(value, i_reduced, j_reduced);
//      }
//      catch (const std::out_of_range& e) {
//      }
//   }
//   return coo_matrix;
//}

//CSCSymmetricMatrix UnoMatrix::to_CSC() {
//   CSCSymmetricMatrix csc_matrix(this->dimension);
//
//   int current_column = 0;
//   int number_terms = 0;
//   csc_matrix.column_start.push_back(0);
//   for (const auto[key, value]: this->matrix) {
//      csc_matrix.matrix.push_back(value);
//      // retrieve indices
//      int i = key % this->dimension;
//      int j = key / this->dimension;
//
//      csc_matrix.row_number.push_back(i);
//      for (int column = current_column; column < j; column++) {
//         csc_matrix.column_start.push_back(number_terms);
//      }
//      current_column = j;
//      number_terms++;
//   }
//   csc_matrix.column_start.push_back(number_terms);
//   return csc_matrix;
//}

//CSCSymmetricMatrix UnoMatrix::to_CSC(const std::unordered_map<int, int>& mask) {
//   CSCSymmetricMatrix csc_matrix(this->dimension);
//
//   int current_column = 0;
//   int number_terms = 0;
//   csc_matrix.column_start.push_back(0);
//   for (const auto[key, value]: this->matrix) {
//      // retrieve indices
//      int i = key % this->dimension;
//      int j = key / this->dimension;
//      try {
//         // if i and j are kept, compute the indices in the reduced matrix
//         int i_reduced = mask.at(i);
//         int j_reduced = mask.at(j);
//         csc_matrix.matrix.push_back(value);
//         csc_matrix.row_number.push_back(i_reduced);
//         for (int column = current_column; column < j_reduced; column++) {
//            csc_matrix.column_start.push_back(number_terms);
//         }
//         current_column = j_reduced;
//         number_terms++;
//      }
//      catch (const std::out_of_range& e) {
//      }
//   }
//   csc_matrix.column_start.push_back(number_terms);
//   return csc_matrix;
//}
//
//std::ostream& operator<<(std::ostream& stream, UnoMatrix& matrix) {
//   for (const auto[key, value]: matrix.matrix) {
//      // retrieve indices
//      size_t i = key % matrix.dimension;
//      size_t j = key / matrix.dimension;
//      stream << "m(" << i << ", " << j << ") = " << value << ", ";
//   }
//   return stream;
//}

