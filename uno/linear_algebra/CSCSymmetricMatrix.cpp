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
// padding provides padding_size slots at the end of each column

CSCSymmetricMatrix::CSCSymmetricMatrix(int dimension, size_t capacity, size_t padding_size) : SymmetricMatrix(dimension, capacity),
matrix(capacity + dimension*padding_size), column_start(dimension + 1), row_index(capacity + dimension*padding_size),
remaining_column_padding(dimension, padding_size) {
}

CSCSymmetricMatrix::CSCSymmetricMatrix(std::vector<double> matrix, const std::vector<int>& column_start, std::vector<int> row_number, int
capacity) : SymmetricMatrix((int) column_start.size() - 1, capacity), matrix(std::move(matrix)), column_start(column_start), row_index(std::move
(row_number)), remaining_column_padding(dimension, 0) {
   assert(false && "padding should be fixed");
}

// generic iterator
void CSCSymmetricMatrix::for_each(const std::function<void (int, int, double)>& f) const {
   size_t overall_padding_size = 0;
   for (int j = 0; j < this->dimension; j++) {
      // work on the active elements (ignore the padding slots)
      for (int k = this->column_start[j] + overall_padding_size; k < this->column_start[j + 1] + overall_padding_size; k++) {
         int i = this->row_index[k];
         f(i, j, this->matrix[k]);
      }
      overall_padding_size += this->remaining_column_padding[j];
   }
}

void CSCSymmetricMatrix::insert(double term, int row_index, int column_index) {
   assert(this->current_insertion_index <= this->matrix.capacity() && "The matrix does not have sufficient capacity");
   assert(column_index == this->current_column && "The previous columns should be finalized");
   this->matrix[this->current_insertion_index] = term;
   this->row_index[this->current_insertion_index] = row_index;
   this->column_start[column_index + 1]++;
   this->current_insertion_index++;
   this->number_nonzeros++;
}

void CSCSymmetricMatrix::finalize(size_t column_index) {
   assert(column_index == this->current_column && "You are not finalizing the current column");
   assert(column_index < this->dimension && "The dimension of the matrix was exceeded");
   // add padding
   this->current_insertion_index += this->remaining_column_padding[column_index];

   if (column_index < this->dimension - 1) {
      this->column_start[column_index + 2] = this->column_start[column_index + 1];
      this->current_column++;
   }
}

void CSCSymmetricMatrix::add_identity_multiple(double factor) {
   size_t overall_padding_size = 0;
   for (int j = 0; j < this->dimension; j++) {
      // check if the last element in the column is diagonal
      const size_t number_elements = this->column_start[j+1] - this->column_start[j];
      bool insert_new_element = false;
      if (0 < number_elements) {
         const size_t last_element_index = this->column_start[j+1] - 1 + overall_padding_size;
         const int i = this->row_index[last_element_index];
         if (i == j) {
            // the diagonal element already exists: simply add factor
            this->matrix[last_element_index] += factor;
         }
         else { // column terminates at a non-diagonal element
            insert_new_element = true;
         }
      }
      else { // empty column
         insert_new_element = true;
      }

      // insert a new element: there must be at least one padding slot
      if (insert_new_element) {
         assert(0 < this->remaining_column_padding[j] && "Padding is not sufficient to insert a new element");
         const size_t last_element_index = this->column_start[j+1] - 1 + overall_padding_size;
         // insert the new diagonal element at the first padding slot
         this->matrix[last_element_index+1] = factor;
         this->row_index[last_element_index+1] = j;
         // shift the next column starts
         for (int k = j + 1; k < this->dimension + 1; k++) {
            this->column_start[k]++;
         }
         this->number_nonzeros++;
         this->remaining_column_padding[j]--;
      }
      overall_padding_size += this->remaining_column_padding[j];
   }
}

void CSCSymmetricMatrix::force_explicit_diagonal_elements() {
   this->add_identity_multiple(0.);
}

COOSymmetricMatrix CSCSymmetricMatrix::to_COO() {
   /*
   COOSymmetricMatrix coo_matrix(this->dimension, this->capacity);
   this->for_each([&](int i, int j, double entry) {
      coo_matrix.insert(entry, i, j);
   });
   return coo_matrix;
    */
   COOSymmetricMatrix coo_matrix(this->dimension, this->matrix, this->column_start, this->row_index);
   //this->for_each([&](int i, int j, double entry) {
   //   coo_matrix.insert(entry, i, j);
   //});
   return coo_matrix;
}

CSCSymmetricMatrix CSCSymmetricMatrix::identity(int dimension) {
   // no padding
   CSCSymmetricMatrix identity(dimension, dimension, 0);
   for (int i = 0; i < dimension; i++) {
      identity.insert(1., i, i);
      identity.finalize(i);
   }
   return identity;
}

std::ostream& operator<<(std::ostream& stream, CSCSymmetricMatrix& matrix) {
   stream << matrix.number_nonzeros << " non zeros\n";
   stream << "W =";
   matrix.for_each([&](size_t /*i*/, size_t /*j*/, double entry) {
      stream << " " << entry;
   });
   stream << "\n";
   stream << "with column start: ";
   print_vector(stream, matrix.column_start, '\n', 0, matrix.dimension + 1);
   stream << "and row index:";
   matrix.for_each([&](size_t i, size_t /*j*/, double /*entry*/) {
      stream << " " << i;
   });
   stream << "\n";
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