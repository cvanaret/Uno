#include <ostream>
#include <cassert>
#include <utility>
#include "CSCSymmetricMatrix.hpp"
#include "Vector.hpp"

/*
 * Compressed Sparse Column
 * https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
 */

// matrix and row_index have nnz elements
// column_start has dimension+1 elements

CSCSymmetricMatrix::CSCSymmetricMatrix(int dimension, size_t capacity) : SymmetricMatrix(dimension, capacity),
      matrix(capacity), column_start(dimension + 1), row_index(capacity) {
}

CSCSymmetricMatrix::CSCSymmetricMatrix(std::vector<double> matrix, const std::vector<int>& column_start, std::vector<int> row_number, int
capacity) : SymmetricMatrix((int) column_start.size() - 1, capacity), matrix(std::move(matrix)), column_start(column_start), row_index(std::move
(row_number)) {
   //assert(false && "CSCSymmetricMatrix::CSCSymmetricMatrix to check");
}

// generic iterator
void CSCSymmetricMatrix::for_each(const std::function<void (int, int, double)>& f) const {
   for (int j = 0; j < this->dimension; j++) {
      for (int k = this->column_start[j]; k < this->column_start[j + 1]; k++) {
         int i = this->row_index[k];
         f(i, j, this->matrix[k]);
      }
   }
}

void CSCSymmetricMatrix::insert(double /*term*/, int /*row_index*/, int /*column_index*/) {
   assert(false && "CSCSymmetricMatrix::insert is not implemented");
}

CSCSymmetricMatrix CSCSymmetricMatrix::add_identity_multiple(double multiple) {
   /* initialize the damped matrix */
   std::vector<double> damped_matrix;
   std::vector<int> damped_column_start;
   std::vector<int> damped_row_number;
   damped_matrix.reserve(this->capacity);
   damped_row_number.reserve(this->capacity);
   damped_column_start.reserve(this->dimension);

   int current_number_nonzeros = 0;
   damped_column_start.push_back(0);

   /* go through the columns */
   for (int j = 0; j < this->dimension; j++) {
      bool diagonal_term_updated = false;

      for (int k = this->column_start[j]; k < this->column_start[j + 1]; k++) {
         /* compute row number */
         int i = this->row_index[k];

         if (i == j) { /* update diagonal term */
            damped_matrix.push_back(this->matrix[k] + multiple);
            diagonal_term_updated = true;
         }
         else if (j < i && !diagonal_term_updated) { /* we passed the diagonal (j, j) */
            damped_matrix.push_back(multiple);
            damped_matrix.push_back(this->matrix[k]);
            damped_row_number.push_back(j); // diagonal term
            current_number_nonzeros++;
            diagonal_term_updated = true;
         }
         else { /* keep off-diagonal term */
            damped_matrix.push_back(this->matrix[k]);
         }
         damped_row_number.push_back(i);
         current_number_nonzeros++;
      }
      /* add diagonal term in column j if not present */
      if (!diagonal_term_updated) {
         damped_matrix.push_back(multiple);
         damped_row_number.push_back(j);
         current_number_nonzeros++;
      }
      damped_column_start.push_back(current_number_nonzeros);
   }
   return {damped_matrix, damped_column_start, damped_row_number, (int) this->capacity};
}

COOSymmetricMatrix CSCSymmetricMatrix::to_COO() {
   COOSymmetricMatrix coo_matrix(this->dimension, this->capacity);
   this->for_each([&](int i, int j, double entry) {
      coo_matrix.insert(entry, i, j);
   });
   return coo_matrix;
}

CSCSymmetricMatrix CSCSymmetricMatrix::identity(int dimension) {
   /* initialize the identity matrix */
   std::vector<double> matrix(dimension);
   std::vector<int> column_start(dimension + 1);
   std::vector<int> row_index(dimension);

   column_start[0] = 0;
   for (int i = 0; i < dimension; i++) {
      matrix[i] = 1.;
      row_index[i] = i;
      column_start[i + 1] = i + 1;
   }
   return {matrix, column_start, row_index, dimension};
}

std::ostream& operator<<(std::ostream& stream, CSCSymmetricMatrix& matrix) {
   stream << matrix.number_nonzeros << " non zeros\n";
   stream << "W = ";
   print_vector(stream, matrix.matrix, '\n', 0, matrix.number_nonzeros);
   stream << "with column start: ";
   print_vector(stream, matrix.column_start, '\n', 0, matrix.dimension + 1);
   stream << "and row index: ";
   print_vector(stream, matrix.row_index, '\n', 0, matrix.number_nonzeros);
   return stream;
}

std::ostream& operator<<(std::ostream& stream, const CSCSymmetricMatrix& matrix) {
   stream << matrix.number_nonzeros << " non zeros\n";
   stream << "W = ";
   print_vector(stream, matrix.matrix, '\n', 0, matrix.number_nonzeros);
   stream << "with column start: ";
   print_vector(stream, matrix.column_start, '\n', 0, matrix.dimension + 1);
   stream << "and row index: ";
   print_vector(stream, matrix.row_index, '\n', 0, matrix.number_nonzeros);
   return stream;
}