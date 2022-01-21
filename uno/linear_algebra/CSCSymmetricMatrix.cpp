#include <ostream>
#include <cassert>
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
column_starts(dimension + 1), remaining_column_padding(dimension, padding_size) {
   const size_t capacity_with_padding = capacity + dimension*padding_size;
   entries.reserve(capacity_with_padding);
   row_indices.reserve(capacity_with_padding);
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
   size_t overall_padding_size = 0;
   for (size_t j = 0; j < this->dimension; j++) {
      // work on the active elements (ignore the padding slots)
      for (size_t k = this->column_starts[j] + overall_padding_size; k < this->column_starts[j + 1] + overall_padding_size; k++) {
         const size_t i = this->row_indices[k];
         const double entry = this->entries[k];
         f(i, j, entry);
      }
      overall_padding_size += this->remaining_column_padding[j];
   }
}

void CSCSymmetricMatrix::for_each(size_t column_index, const std::function<void (size_t, double)>& f) const {
   size_t initial_value = 0;
   const size_t overall_padding_size = std::accumulate(std::begin(this->remaining_column_padding), std::next(std::begin
      (this->remaining_column_padding), column_index), initial_value);
   // work on the active elements (ignore the padding slots)
   for (size_t k = this->column_starts[column_index] + overall_padding_size; k < this->column_starts[column_index + 1] + overall_padding_size; k++) {
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

void CSCSymmetricMatrix::pop() {
   assert(false && "CSC::pop not implemented");
}

void CSCSymmetricMatrix::finalize(size_t column_index) {
   assert(column_index == this->current_column && "You are not finalizing the current column");
   assert(column_index < this->dimension && "The dimension of the matrix was exceeded");

   // add padding
   for (size_t k = 0; k < this->remaining_column_padding[column_index]; k++) {
      this->entries.push_back(0.);
      this->row_indices.push_back(0);
   }

   // start the next column at the current start
   if (column_index < this->dimension - 1) {
      this->column_starts[column_index + 2] = this->column_starts[column_index + 1];
      this->current_column++;
   }
}

void CSCSymmetricMatrix::add_identity_multiple(double multiple) {
   size_t overall_padding_size = 0;
   // go through each column
   for (size_t j = 0; j < this->dimension; j++) {
      // check if the last element in the column is diagonal
      const size_t number_elements = this->column_starts[j + 1] - this->column_starts[j];
      // a new element will be inserted unless it already exists
      bool insert_new_element = true;
      if (0 < number_elements) {
         const size_t last_entry_index = this->column_starts[j + 1] - 1 + overall_padding_size;
         const size_t i = this->row_indices[last_entry_index];
         if (i == j) {
            // the diagonal element already exists: simply add factor
            this->entries[last_entry_index] += multiple;
            insert_new_element = false;
         }
         // otherwise, column terminates at a non-diagonal element
      }

      // insert a new element: there must be at least one padding slot
      if (insert_new_element) {
         assert(0 < this->remaining_column_padding[j] && "Padding is not sufficient to insert a new element");
         const size_t last_element_index = this->column_starts[j + 1] - 1 + overall_padding_size;
         // insert the new diagonal element at the first padding slot
         this->entries[last_element_index + 1] = multiple;
         this->row_indices[last_element_index + 1] = j;
         // shift the next column starts
         for (size_t k = j + 1; k < this->dimension + 1; k++) {
            this->column_starts[k]++;
         }
         this->number_nonzeros++;
         this->remaining_column_padding[j]--;
      }
      overall_padding_size += this->remaining_column_padding[j];
   }
}

double CSCSymmetricMatrix::smallest_diagonal_entry() const {
   double smallest_entry = std::numeric_limits<double>::infinity();
   size_t overall_padding_size = 0;
   // go through each column
   for (size_t j = 0; j < this->dimension; j++) {
      // if it exists, the diagonal entry is the last element of the column
      const size_t number_elements = this->column_starts[j + 1] - this->column_starts[j];
      if (0 < number_elements) {
         const size_t last_entry_index = this->column_starts[j + 1] - 1 + overall_padding_size;
         // get the row index of the last column entry
         const size_t i = this->row_indices[last_entry_index];
         // if the entry is diagonal, update the smallest diagonal entry
         if (i == j) {
            smallest_entry = std::min(smallest_entry, this->entries[last_entry_index]);
         }
      }
      overall_padding_size += this->remaining_column_padding[j];
   }
   if (smallest_entry == std::numeric_limits<double>::infinity()) {
      smallest_entry = 0.;
   }
   return smallest_entry;
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

CSCSymmetricMatrix CSCSymmetricMatrix::identity(size_t dimension) {
   // no padding
   CSCSymmetricMatrix identity(dimension, dimension, 0);
   for (size_t i = 0; i < dimension; i++) {
      identity.insert(1., i, i);
      identity.finalize(i);
   }
   return identity;
}

void CSCSymmetricMatrix::print(std::ostream& stream) const {
   stream << "W = ";
   print_vector(stream, this->entries, 0, this->number_nonzeros);
   stream << "with column start: ";
   print_vector(stream, this->column_starts, 0, this->dimension + 1);
   stream << "and row index: ";
   print_vector(stream, this->row_indices, 0, this->number_nonzeros);
}