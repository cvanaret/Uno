#include <ostream>
#include <cassert>
#include <utility>
#include <numeric>
#include "CSCSymmetricMatrix.hpp"
#include "Vector.hpp"

/*
 * Compressed Sparse Column
 * https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
 */

// matrix and row_index have nnz elements
// column_start has dimension+1 elements
// padding provides padding_size slots at the end of each column

CSCSymmetricMatrix::CSCSymmetricMatrix(size_t dimension, size_t capacity, size_t padding_size) : SymmetricMatrix(dimension, capacity),
column_start(dimension + 1), remaining_column_padding(dimension, padding_size) {
   const size_t capacity_with_padding = capacity + dimension*padding_size;
   matrix.reserve(capacity_with_padding);
   row_index.reserve(capacity_with_padding);
}

CSCSymmetricMatrix::CSCSymmetricMatrix(std::vector<double> matrix, const std::vector<int>& column_start, std::vector<int> row_number, int
capacity) : SymmetricMatrix(column_start.size() - 1, capacity), matrix(std::move(matrix)), column_start(column_start), row_index(std::move
(row_number)), remaining_column_padding(dimension, 0) {
   assert(false && "padding should be fixed");
}

// generic iterator
void CSCSymmetricMatrix::for_each(const std::function<void (int, int, double)>& f) const {
   size_t overall_padding_size = 0;
   for (size_t j = 0; j < this->dimension; j++) {
      // work on the active elements (ignore the padding slots)
      for (size_t k = this->column_start[j] + overall_padding_size; k < this->column_start[j + 1] + overall_padding_size; k++) {
         int i = this->row_index[k];
         f(i, static_cast<int>(j), this->matrix[k]);
      }
      overall_padding_size += this->remaining_column_padding[j];
   }
}

void CSCSymmetricMatrix::for_each(size_t column_index, const std::function<void (int, double)>& f) const {
   const size_t overall_padding_size = std::accumulate(begin(this->remaining_column_padding), begin(this->remaining_column_padding) + column_index,
         0);
   // work on the active elements (ignore the padding slots)
   for (size_t k = this->column_start[column_index] + overall_padding_size; k < this->column_start[column_index + 1] + overall_padding_size; k++) {
      int i = this->row_index[k];
      f(i, this->matrix[k]);
   }
}

void CSCSymmetricMatrix::insert(double term, size_t row_index, size_t column_index) {
   assert(column_index == this->current_column && "The previous columns should be finalized");

   this->matrix.push_back(term);
   this->row_index.push_back(static_cast<int>(row_index));
   this->column_start[column_index + 1]++;
   this->number_nonzeros++;
}

void CSCSymmetricMatrix::finalize(size_t column_index) {
   assert(column_index == this->current_column && "You are not finalizing the current column");
   assert(column_index < this->dimension && "The dimension of the matrix was exceeded");
   // add padding
   for (size_t k = 0; k < this->remaining_column_padding[column_index]; k++) {
      this->matrix.push_back(0);
      this->row_index.push_back(-1);
   }

   // start the next column at the current start
   if (column_index < this->dimension - 1) {
      this->column_start[column_index + 2] = this->column_start[column_index + 1];
      this->current_column++;
   }
}

void CSCSymmetricMatrix::add_identity_multiple(double factor) {
   size_t overall_padding_size = 0;
   // go through each column
   for (size_t j = 0; j < this->dimension; j++) {
      // check if the last element in the column is diagonal
      const size_t number_elements = this->column_start[j+1] - this->column_start[j];
      bool insert_new_element = false;
      if (0 < number_elements) {
         const size_t last_element_index = this->column_start[j+1] - 1 + overall_padding_size;
         const int i = this->row_index[last_element_index];
         if (i == static_cast<int>(j)) {
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
         this->row_index[last_element_index+1] = static_cast<int>(j);
         // shift the next column starts
         for (size_t k = j + 1; k < this->dimension + 1; k++) {
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

void CSCSymmetricMatrix::remove_variables(const std::vector<int>& /*variable_indices*/) {
   assert(false && "Not yet implemented");
   /*
   assert(in_increasing_order(variable_indices.data(), variable_indices.size()) && "The list of variables is not sorted");

   // remove columns
   for (int j: variable_indices) {
      const size_t number_elements = this->column_start[j+1] - this->column_start[j];
      // move all the terms to the padding
      this->remaining_column_padding[j] += number_elements;
      this->number_nonzeros -= number_elements;
      print_vector(std::cout, this->column_start);
      for (size_t k = j; k < this->dimension; k++) {
         this->column_start[k] = this->column_start[k+1] - number_elements;
      }
      print_vector(std::cout, this->column_start);
      this->dimension--;
   }

   // remove rows
    */
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

CSCSymmetricMatrix CSCSymmetricMatrix::identity(size_t dimension) {
   // no padding
   CSCSymmetricMatrix identity(dimension, dimension, 0);
   for (size_t i = 0; i < dimension; i++) {
      identity.insert(1., i, i);
      identity.finalize(i);
   }
   return identity;
}

std::ostream& operator<<(std::ostream& stream, const CSCSymmetricMatrix& matrix) {
   stream << matrix.dimension << " variables\n";
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