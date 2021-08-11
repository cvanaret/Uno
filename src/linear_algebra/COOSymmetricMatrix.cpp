#include <ostream>
#include "COOSymmetricMatrix.hpp"

/*
 * Coordinate list
 * https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)
 */

COOSymmetricMatrix::COOSymmetricMatrix(int dimension, size_t capacity) : SymmetricMatrix(dimension, capacity) {
   this->matrix.reserve(capacity);
   this->row_indices.reserve(capacity);
   this->column_indices.reserve(capacity);
}

// generic iterator
void COOSymmetricMatrix::for_each(const std::function<void (int, int, double)>& f) const {
   for (int k = 0; k < this->number_nonzeros; k++) {
      int i = this->row_indices[k];
      int j = this->column_indices[k];
      f(i, j, this->matrix[k]);
   }
}

void COOSymmetricMatrix::insert(double term, int row_index, int column_index) {
   //assert(this->number_nonzeros <= this->capacity && "The capacity of the matrix is too low for insertion");

   this->matrix.push_back(term);
   this->row_indices.push_back((int) row_index);
   this->column_indices.push_back((int) column_index);
   this->number_nonzeros++;
}

double COOSymmetricMatrix::norm_1() {
   // compute maximum column index
   int number_columns = 0;
   for (int j : this->column_indices) {
      number_columns = std::max(number_columns, 1 + j);
   }
   // read the matrix and fill in the column_vectors norm vector
   std::vector<double> column_vectors(number_columns);
   for (size_t k = 0; k < this->matrix.size(); k++) {
      int j = this->column_indices[k];
      column_vectors[j] += std::abs(this->matrix[k]);
   }
   // compute the maximal component of the column_vectors vector
   double norm = 0.;
   for (double j : column_vectors) {
      norm = std::max(norm, j);
   }
   return norm;
}

std::ostream& operator<<(std::ostream& stream, COOSymmetricMatrix& matrix) {
   matrix.for_each([&](int i, int j, double entry) {
      stream << "m(" << i << ", " << j << ") = " << entry << "\n";
   });
   return stream;
}

std::ostream& operator<<(std::ostream& stream, const COOSymmetricMatrix& matrix) {
   matrix.for_each([&](int i, int j, double entry) {
      stream << "m(" << i << ", " << j << ") = " << entry << "\n";
   });
   return stream;
}