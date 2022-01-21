#include <exception>
#include <cassert>
#include <iomanip>
#include "SymmetricMatrix.hpp"

SymmetricMatrix::SymmetricMatrix(size_t dimension, size_t capacity) : dimension(dimension), capacity(capacity) {
   entries.reserve(capacity);
}

void SymmetricMatrix::reset() {
   this->number_nonzeros = 0;
}

std::vector<double> SymmetricMatrix::product(const std::vector<double>& vector) const {
   assert(this->dimension == vector.size() && "The matrix and the vector do not have the same size");

   // create (n, 1) result
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

double SymmetricMatrix::quadratic_product(const std::vector<double>& x, const std::vector<double>& y, size_t block_size) const {
   assert(x.size() == y.size() && "SymmetricMatrix::quadratic_product: the two vectors x and y do not have the same size");
   assert(block_size <= x.size() && "SymmetricMatrix::quadratic_product: the block size is larger than the vectors");

   double result = 0.;
   this->for_each([&](size_t i, size_t j, double entry) {
      if (i < block_size && j < block_size) {
         result += (i == j ? 1 : 2) * entry * x[i] * y[j];
      }
   });
   return result;
}

void SymmetricMatrix::add_outer_product(const SparseVector<double>& x, double scaling_factor) {
   // perform matrix addition: A + rho x x^T
   x.for_each([&](size_t row_index, double row_entry) {
      x.for_each([&](size_t column_index, double column_entry) {
         // upper triangular matrix
         if (row_index <= column_index) {
            // add product of components
            this->insert(scaling_factor * row_entry * column_entry, row_index, column_index);
         }
      });
   });
}

void SymmetricMatrix::finalize(size_t /*column_index*/) {
   // by default, do nothing
}

std::ostream& operator<<(std::ostream& stream, const SymmetricMatrix& matrix) {
   stream << matrix.dimension << " variables, " << matrix.number_nonzeros << " non zeros:\n";
   matrix.print(stream);
   return stream;
}