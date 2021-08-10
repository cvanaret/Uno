#include <ostream>
#include <cassert>
#include "CSRSymmetricMatrix.hpp"
#include "Vector.hpp"

/*
 * Compressed Sparse Row
 * https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format)
 */

// matrix and column_index have nnz elements
// row_start has dimension+1 elements

CSRSymmetricMatrix::CSRSymmetricMatrix(size_t dimension, size_t capacity) : SymmetricMatrix(dimension, capacity),
      matrix(capacity), row_start(dimension + 1), column_index(capacity) {
}

CSRSymmetricMatrix::CSRSymmetricMatrix(const std::vector<double>& matrix, const std::vector<size_t>& row_start, const std::vector<size_t>&
column_index, size_t capacity) : SymmetricMatrix(row_start.size() - 1, capacity), matrix(matrix), row_start(row_start), column_index
      (column_index) {
   assert(false && "CSRSymmetricMatrix::CSRSymmetricMatrix to check");
}

// generic iterator
void CSRSymmetricMatrix::for_each(const std::function<void (size_t, size_t, double)>& f) const {
   for (size_t i = 0; i < this->dimension; i++) {
      for (size_t k = this->row_start[i]; k < this->row_start[i + 1]; k++) {
         size_t j = this->column_index[k];
         f(i, j, this->matrix[k]);
      }
   }
}

void CSRSymmetricMatrix::insert(double /*term*/, size_t /*row_index*/, size_t /*column_index*/) {
   assert(false && "CSRSymmetricMatrix::insert is not implemented");
}

CSRSymmetricMatrix CSRSymmetricMatrix::add_identity_multiple(double /*multiple*/) {
   assert(false && "CSRSymmetricMatrix::add_identity_multiple not implemented");
   //   /* initialize the damped matrix */
   //   std::vector<double> damped_matrix;
   //   std::vector<size_t> damped_column_start;
   //   std::vector<size_t> damped_row_number;
   //   damped_matrix.reserve(this->capacity);
   //   damped_row_number.reserve(this->capacity);
   //   damped_column_start.reserve(this->dimension);
   //
   //   int current_number_nonzeros = 0;
   //   damped_column_start.push_back(0);
   //
   //   /* go through the columns */
   //   for (size_t j = 0; j < this->dimension; j++) {
   //      bool diagonal_term_updated = false;
   //
   //      for (size_t k = this->column_start[j]; k < this->column_start[j + 1]; k++) {
   //         /* compute row number */
   //         size_t i = this->row_index[k];
   //
   //         if (i == j) { /* update diagonal term */
   //            damped_matrix.push_back(this->matrix[k] + multiple);
   //            diagonal_term_updated = true;
   //         }
   //         else if (j < i && !diagonal_term_updated) { /* we passed the diagonal (j, j) */
   //            damped_matrix.push_back(multiple);
   //            damped_matrix.push_back(this->matrix[k]);
   //            damped_row_number.push_back(j); // diagonal term
   //            current_number_nonzeros++;
   //            diagonal_term_updated = true;
   //         }
   //         else { /* keep off-diagonal term */
   //            damped_matrix.push_back(this->matrix[k]);
   //         }
   //         damped_row_number.push_back(i);
   //         current_number_nonzeros++;
   //      }
   //      /* add diagonal term in column j if not present */
   //      if (!diagonal_term_updated) {
   //         damped_matrix.push_back(multiple);
   //         damped_row_number.push_back(j);
   //         current_number_nonzeros++;
   //      }
   //      damped_column_start.push_back(current_number_nonzeros);
   //   }
   //   return CSRSymmetricMatrix(damped_matrix, damped_column_start, damped_row_number, this->capacity);
}

COOSymmetricMatrix CSRSymmetricMatrix::to_COO() {
   COOSymmetricMatrix coo_matrix(this->dimension, this->capacity);
   this->for_each([&](size_t i, size_t j, double entry) {
      coo_matrix.insert(entry, i, j);
   });
   return coo_matrix;
}

CSRSymmetricMatrix CSRSymmetricMatrix::identity(size_t dimension) {
   /* initialize the identity matrix */
   std::vector<double> matrix(dimension);
   std::vector<size_t> row_start(dimension + 1);
   std::vector<size_t> column_index(dimension);

   row_start[0] = 0;
   for (size_t i = 0; i < dimension; i++) {
      matrix[i] = 1.;
      row_start[i + 1] = i + 1;
      column_index[i] = i;
   }

   return CSRSymmetricMatrix(matrix, row_start, column_index, dimension);
}

std::ostream& operator<<(std::ostream& stream, CSRSymmetricMatrix& matrix) {
   stream << matrix.number_nonzeros << " non zeros\n";
   stream << "W = ";
   print_vector(stream, matrix.matrix, '\n', 0, matrix.number_nonzeros);
   stream << "with row start: ";
   print_vector(stream, matrix.row_start, '\n', 0, matrix.dimension + 1);
   stream << "and column index: ";
   print_vector(stream, matrix.column_index, '\n', 0, matrix.number_nonzeros);
   return stream;
}

std::ostream& operator<<(std::ostream& stream, const CSRSymmetricMatrix& matrix) {
   stream << matrix.number_nonzeros << " non zeros\n";
   stream << "W = ";
   print_vector(stream, matrix.matrix, '\n', 0, matrix.number_nonzeros);
   stream << "with row start: ";
   print_vector(stream, matrix.row_start, '\n', 0, matrix.dimension + 1);
   stream << "and column index: ";
   print_vector(stream, matrix.column_index, '\n', 0, matrix.number_nonzeros);
   return stream;
}