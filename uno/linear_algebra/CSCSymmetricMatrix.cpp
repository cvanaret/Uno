// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <ostream>
#include <cassert>
#include "CSCSymmetricMatrix.hpp"
#include "Vector.hpp"
#include "tools/Infinity.hpp"

/*
 * Compressed Sparse Column
 * https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
 */

// matrix and row_index have nnz elements
// column_start has dimension+1 elements

CSCSymmetricMatrix::CSCSymmetricMatrix(size_t dimension, size_t original_capacity, bool use_regularization):
      SymmetricMatrix(dimension, original_capacity, use_regularization),
      column_starts(dimension + 1) {
   entries.reserve(this->capacity);
   row_indices.reserve(this->capacity);
}

void CSCSymmetricMatrix::reset() {
   SymmetricMatrix::reset();
   this->entries.clear();
   this->row_indices.clear();
   initialize_vector<size_t>(this->column_starts, 0);
   this->current_column = 0;
}

// generic iterator
void CSCSymmetricMatrix::for_each(const std::function<void (size_t, size_t, double)>& f) const {
   for (size_t j = 0; j < this->dimension; j++) {
      for (size_t k = this->column_starts[j]; k < this->column_starts[j + 1]; k++) {
         const size_t i = this->row_indices[k];
         const double entry = this->entries[k];
         f(i, j, entry);
      }
   }
}

void CSCSymmetricMatrix::for_each(size_t column_index, const std::function<void (size_t, double)>& f) const {
   for (size_t k = this->column_starts[column_index]; k < this->column_starts[column_index + 1]; k++) {
      const size_t i = this->row_indices[k];
      const double entry = this->entries[k];
      f(i, entry);
   }
}

void CSCSymmetricMatrix::insert(double term, size_t row_index, size_t column_index) {
   assert(column_index == this->current_column && "The previous columns should be finalized");

   this->entries.push_back(term);
   this->row_indices.push_back(row_index);
   this->column_starts[column_index + 1]++;
   this->number_nonzeros++;
}

void CSCSymmetricMatrix::finalize_column(size_t column_index) {
   assert(column_index == this->current_column && "You are not finalizing the current column");
   assert(column_index < this->dimension && "The dimension of the matrix was exceeded");

   // possibly add regularization
   if (this->use_regularization) {
      insert(0., column_index, column_index);
   }
   this->current_column++;

   // start the next column at the current start
   if (column_index < this->dimension - 1) {
      this->column_starts[column_index + 2] = this->column_starts[column_index + 1];
   }
}

double CSCSymmetricMatrix::smallest_diagonal_entry() const {
   double smallest_entry = INF;
   // go through each column
   for (size_t j = 0; j < this->dimension; j++) {
      // if it exists, the diagonal entry is the last element of the column
      const size_t number_elements = this->column_starts[j + 1] - this->column_starts[j];
      if (0 < number_elements) {
         const size_t last_entry_index = this->column_starts[j + 1] - 1;
         // get the row index of the last column entry
         const size_t i = this->row_indices[last_entry_index];
         // if the entry is diagonal, update the smallest diagonal entry
         if (i == j) {
            smallest_entry = std::min(smallest_entry, this->entries[last_entry_index]);
         }
      }
   }
   if (smallest_entry == INF) {
      smallest_entry = 0.;
   }
   return smallest_entry;
}

void CSCSymmetricMatrix::set_regularization(const std::function<double(size_t index)>& regularization_function) {
   assert(this->use_regularization && "You are trying to regularize a matrix where regularization was not preallocated.");
   for (size_t i = 0; i < this->dimension; i++) {
      // the regularization term is located at the end of the column, that is right before the start of the next column
      const size_t k = this->column_starts[i + 1] - 1;
      this->entries[k] = regularization_function(i);
   }
}

void CSCSymmetricMatrix::print(std::ostream& stream) const {
   stream << "W = ";
   print_vector(stream, this->entries, 0, this->number_nonzeros);
   stream << "with column start: ";
   print_vector(stream, this->column_starts, 0, this->dimension + 1);
   stream << "and row index: ";
   print_vector(stream, this->row_indices, 0, this->number_nonzeros);
}