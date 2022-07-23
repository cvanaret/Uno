// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <iostream>
#include "COOSymmetricMatrix.hpp"
#include "tools/Infinity.hpp"

/*
 * Coordinate list
 * https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)
 */

COOSymmetricMatrix::COOSymmetricMatrix(size_t dimension, size_t original_capacity, bool use_regularization):
      SymmetricMatrix(dimension, original_capacity, use_regularization) {
   row_indices.reserve(this->capacity);
   column_indices.reserve(this->capacity);

   if (this->use_regularization) {
      this->initialize_regularization();
   }
}

void COOSymmetricMatrix::reset() {
   // empty the matrix
   SymmetricMatrix::reset();
   this->entries.clear();
   this->row_indices.clear();
   this->column_indices.clear();

   // initialize regularization terms
   if (this->use_regularization) {
      this->initialize_regularization();
   }
}

// generic iterator
void COOSymmetricMatrix::for_each(const std::function<void(size_t, size_t, double)>& f) const {
   for (size_t k = 0; k < this->number_nonzeros; k++) {
      const size_t i = this->row_indices[k];
      const size_t j = this->column_indices[k];
      const double entry = this->entries[k];
      f(i, j, entry);
   }
}

void COOSymmetricMatrix::insert(double term, size_t row_index, size_t column_index) {
   assert(this->number_nonzeros <= row_indices.size() && "The COO matrix doesn't have a sufficient capacity");
   this->entries.push_back(term);
   this->row_indices.push_back(row_index);
   this->column_indices.push_back(column_index);
   this->number_nonzeros++;
}

void COOSymmetricMatrix::finalize_column(size_t /*column_index*/) {
   // do nothing
}

double COOSymmetricMatrix::smallest_diagonal_entry() const {
   // TODO: there may be several entries for given indices
   double smallest_entry = INF;
   this->for_each([&](size_t i, size_t j, double entry) {
      if (i == j) {
         smallest_entry = std::min(smallest_entry, entry);
      }
   });
   if (smallest_entry == INF) {
      smallest_entry = 0.;
   }
   return smallest_entry;
}

void COOSymmetricMatrix::set_regularization(const std::function<double(size_t index)>& regularization_function) {
   assert(this->use_regularization && "You are trying to regularize a matrix where regularization was not preallocated.");
   // the regularization terms (that lie at the start of the entries vector) can be directly modified
   for (size_t i = 0; i < this->dimension; i++) {
      this->entries[i] = regularization_function(i);
   }
}

void COOSymmetricMatrix::print(std::ostream& stream) const {
   this->for_each([&](size_t i, size_t j, double entry) {
      stream << "m(" << i << ", " << j << ") = " << entry << '\n';
   });
}

void COOSymmetricMatrix::initialize_regularization() {
   // introduce elements at the start of the entries
   for (size_t i = 0; i < this->dimension; i++) {
      this->insert(0., i, i);
   }
}