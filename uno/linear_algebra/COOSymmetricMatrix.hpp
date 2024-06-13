// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_COOSYMMETRICMATRIX_H
#define UNO_COOSYMMETRICMATRIX_H

#include <ostream>
#include <cassert>
#include "SymmetricMatrix.hpp"
#include "tools/Infinity.hpp"

/*
 * Coordinate list
 * https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)
 */
template <typename ElementType>
class COOSymmetricMatrix: public SymmetricMatrix<ElementType> {
public:
   COOSymmetricMatrix(size_t dimension, size_t original_capacity, bool use_regularization);

   void reset() override;
   void insert(ElementType term, size_t row_index, size_t column_index) override;
   void finalize_column(size_t /*column_index*/) override { /* do nothing */ }
   [[nodiscard]] ElementType smallest_diagonal_entry(size_t max_dimension) const override;
   void set_regularization(const std::function<ElementType(size_t index)>& regularization_function) override;

   void print(std::ostream& stream) const override;

   static COOSymmetricMatrix<ElementType> zero(size_t dimension);

protected:
   std::vector<size_t> row_indices;
   std::vector<size_t> column_indices;
   Vector<ElementType> diagonal_entries;

   void initialize_regularization();

   // iterator functions
   [[nodiscard]] std::tuple<size_t, size_t, ElementType> dereference_iterator(size_t column_index, size_t nonzero_index) const override;
   void increment_iterator(size_t& column_index, size_t& nonzero_index) const override;
};

// implementation

template <typename ElementType>
COOSymmetricMatrix<ElementType>::COOSymmetricMatrix(size_t dimension, size_t original_capacity, bool use_regularization):
      SymmetricMatrix<ElementType>(dimension, original_capacity, use_regularization),
      diagonal_entries(dimension, ElementType(0)) {
   this->row_indices.reserve(this->capacity);
   this->column_indices.reserve(this->capacity);

   // initialize regularization terms
   if (this->use_regularization) {
      this->initialize_regularization();
   }
}

template <typename ElementType>
void COOSymmetricMatrix<ElementType>::reset() {
   // empty the matrix
   SymmetricMatrix<ElementType>::reset();
   this->row_indices.clear();
   this->column_indices.clear();
   this->diagonal_entries.fill(ElementType(0));

   // initialize regularization terms
   if (this->use_regularization) {
      this->initialize_regularization();
   }
}

template <typename ElementType>
void COOSymmetricMatrix<ElementType>::insert(ElementType term, size_t row_index, size_t column_index) {
   assert(this->number_nonzeros <= this->row_indices.size() && "The COO matrix doesn't have a sufficient capacity");
   this->entries.push_back(term);
   this->row_indices.push_back(row_index);
   this->column_indices.push_back(column_index);
   this->number_nonzeros++;

   // possibly update diagonal
   if (row_index == column_index) {
      this->diagonal_entries[row_index] += term;
   }
}

template <typename ElementType>
ElementType COOSymmetricMatrix<ElementType>::smallest_diagonal_entry(size_t max_dimension) const {
   ElementType smallest_entry = INF<ElementType>;
   for (size_t row_index: Range(std::min(this->dimension, max_dimension))) {
      smallest_entry = std::min(smallest_entry, this->diagonal_entries[row_index]);
   }
   return smallest_entry;
}

template <typename ElementType>
void COOSymmetricMatrix<ElementType>::set_regularization(const std::function<ElementType(size_t /*index*/)>& regularization_function) {
   assert(this->use_regularization && "You are trying to regularize a matrix where regularization was not preallocated.");

   // the regularization terms (that lie at the start of the entries vector) can be directly modified
   for (size_t row_index: Range(this->dimension)) {
      const ElementType element = regularization_function(row_index);
      this->entries[row_index] = element;
      // update diagonal
      this->diagonal_entries[row_index] += element;
   }
}

template <typename ElementType>
void COOSymmetricMatrix<ElementType>::print(std::ostream& stream) const {
   for (const auto [row_index, column_index, element]: *this) {
      stream << "m(" << row_index << ", " << column_index << ") = " << element << '\n';
   }
}

template <typename ElementType>
void COOSymmetricMatrix<ElementType>::initialize_regularization() {
   // introduce elements at the start of the entries
   for (size_t row_index: Range(this->dimension)) {
      this->insert(ElementType(0), row_index, row_index);
   }
}

template <typename ElementType>
std::tuple<size_t, size_t, ElementType> COOSymmetricMatrix<ElementType>::dereference_iterator(size_t /*column_index*/, size_t nonzero_index) const {
   return {this->row_indices[nonzero_index], this->column_indices[nonzero_index], this->entries[nonzero_index]};
}

template <typename ElementType>
void COOSymmetricMatrix<ElementType>::increment_iterator(size_t& column_index, size_t& nonzero_index) const {
   nonzero_index++;
   // if end reached
   if (nonzero_index == this->number_nonzeros) {
      column_index = this->dimension;
   }
}

template <typename ElementType>
COOSymmetricMatrix<ElementType> COOSymmetricMatrix<ElementType>::zero(size_t dimension) {
   return COOSymmetricMatrix<ElementType>(dimension, 0, false);
}

#endif // UNO_COOSYMMETRICMATRIX_H