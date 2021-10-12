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

COOSymmetricMatrix::COOSymmetricMatrix(size_t dimension, std::vector<double> matrix, std::vector<size_t> row_indices, std::vector<size_t>
      column_indices) :
      SymmetricMatrix(dimension, matrix.capacity()), matrix(std::move(matrix)), row_indices(std::move(row_indices)),
      column_indices(std::move(column_indices)) {
   assert(false && "COOSymmetricMatrix constructor 2");
}

// generic iterator
void COOSymmetricMatrix::for_each(const std::function<void(size_t, size_t, double)>& f) const {
   for (size_t k = 0; k < this->number_nonzeros; k++) {
      size_t i = this->row_indices[k];
      size_t j = this->column_indices[k];
      f(i, j, this->matrix[k]);
   }
}

size_t COOSymmetricMatrix::find(size_t row_index, size_t column_index) {
   size_t position = 0;
   while (position < this->number_nonzeros) {
      if (this->row_indices[position] == row_index && this->column_indices[position] == column_index) {
         break;
      }
      position++;
   }
   return position;
}

void COOSymmetricMatrix::insert(double term, size_t row_index, size_t column_index) {
   /*
   size_t position = this->find(row_index, column_index);
   if (position < this->number_nonzeros) {
      this->matrix[position] += term;
   }
   else {
    */
      this->matrix[this->number_nonzeros] = term;
      this->row_indices[this->number_nonzeros] = row_index;
      this->column_indices[this->number_nonzeros] = column_index;
      this->number_nonzeros++;
   //}
}

std::ostream& operator<<(std::ostream& stream, const COOSymmetricMatrix& matrix) {
   matrix.for_each([&](size_t i, size_t j, double entry) {
      stream << "m(" << i << ", " << j << ") = " << entry << "\n";
   });
   return stream;
}
void COOSymmetricMatrix::add_identity_multiple(double /*multiple*/) {
   // do nothing
   assert(false && "not yet implemented");
}
