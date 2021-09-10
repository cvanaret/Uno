#include <iostream>
#include <cassert>
#include "COOSymmetricMatrix.hpp"

/*
 * Coordinate list
 * https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)
 */

COOSymmetricMatrix::COOSymmetricMatrix(size_t dimension, size_t capacity) : SymmetricMatrix(dimension, capacity), matrix(capacity),
row_indices(capacity), column_indices(capacity) {
}

COOSymmetricMatrix::COOSymmetricMatrix(size_t dimension, std::vector<double> matrix, std::vector<int> row_indices, std::vector<int> column_indices) :
      SymmetricMatrix(dimension, matrix.capacity()), matrix(std::move(matrix)), row_indices(std::move(row_indices)),
      column_indices(std::move(column_indices)) {
   assert(false && "COOSymmetricMatrix constructor 2");
}

// generic iterator
void COOSymmetricMatrix::for_each(const std::function<void(int, int, double)>& f) const {
   for (size_t k = 0; k < this->number_nonzeros; k++) {
      int i = this->row_indices[k];
      int j = this->column_indices[k];
      f(i, j, this->matrix[k]);
   }
}

void COOSymmetricMatrix::insert(double term, size_t row_index, size_t column_index) {
   this->matrix[this->number_nonzeros] = term;
   this->row_indices[this->number_nonzeros] = static_cast<int>(row_index);
   this->column_indices[this->number_nonzeros] = static_cast<int>(column_index);
   this->number_nonzeros++;
}

std::ostream& operator<<(std::ostream& stream, const COOSymmetricMatrix& matrix) {
   matrix.for_each([&](int i, int j, double entry) {
      stream << "m(" << i << ", " << j << ") = " << entry << "\n";
   });
   return stream;
}
void COOSymmetricMatrix::add_identity_multiple(double /*multiple*/) {
   // do nothing
   assert(false && "not yet implemented");
}
