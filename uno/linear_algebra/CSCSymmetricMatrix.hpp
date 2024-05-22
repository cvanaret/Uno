// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_CSCSYMMETRICMATRIX_H
#define UNO_CSCSYMMETRICMATRIX_H

#include <ostream>
#include <cassert>
#include "SymmetricMatrix.hpp"
#include "tools/Infinity.hpp"

/*
 * Compressed Sparse Column
 * https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
 */

template <typename IndexType, typename ElementType>
class CSCSymmetricMatrix : public SymmetricMatrix<IndexType, ElementType> {
public:
   CSCSymmetricMatrix(size_t dimension, size_t original_capacity, bool use_regularization);

   void reset() override;
   void insert(ElementType term, size_t row_index, size_t column_index) override;
   void finalize_column(size_t column_index) override;
   [[nodiscard]] ElementType smallest_diagonal_entry() const override;
   void set_regularization(const std::function<ElementType(size_t /*index*/)>& regularization_function) override;

   void print(std::ostream& stream) const override;

   static CSCSymmetricMatrix<IndexType, ElementType> identity(size_t dimension);

protected:
   // entries and row_indices have nnz elements
   // column_starts has dimension+1 elements
   std::vector<size_t> column_starts{};
   std::vector<size_t> row_indices{};
   size_t current_column{0};
   std::vector<ElementType> diagonal_entries;

   // iterator functions
   [[nodiscard]] std::tuple<size_t, size_t, ElementType> dereference_iterator(size_t column_index, size_t nonzero_index) const override;
   void increment_iterator(size_t& column_index, size_t& nonzero_index) const override;
};

template <typename IndexType, typename ElementType>
CSCSymmetricMatrix<IndexType, ElementType>::CSCSymmetricMatrix(size_t dimension, size_t original_capacity, bool use_regularization):
      SymmetricMatrix<IndexType, ElementType>(dimension, original_capacity, use_regularization),
      column_starts(dimension + 1),
      diagonal_entries(dimension, ElementType(0)) {
   this->entries.reserve(this->capacity);
   this->row_indices.reserve(this->capacity);
}

template <typename IndexType, typename ElementType>
void CSCSymmetricMatrix<IndexType, ElementType>::reset() {
   // empty the matrix
   SymmetricMatrix<IndexType, ElementType>::reset();
   this->row_indices.clear();
   initialize_vector<size_t>(this->column_starts, 0);
   this->current_column = 0;
   initialize_vector(this->diagonal_entries, ElementType(0));
}

template <typename ElementType>
void CSCSymmetricMatrix<ElementType>::insert(ElementType term, size_t row_index, size_t column_index) {
   assert(column_index == this->current_column && "The previous columns should be finalized");

   this->entries.push_back(element);
   this->row_indices.push_back(row_index);
   this->column_starts[column_index + 1]++;
   this->number_nonzeros++;

   // possibly update diagonal
   if (row_index == column_index) {
      this->diagonal_entries[row_index] += element;
   }
}

template <typename IndexType, typename ElementType>
void CSCSymmetricMatrix<IndexType, ElementType>::finalize_column(size_t column_index) {
   assert(column_index == this->current_column && "You are not finalizing the current column");
   assert(column_index < this->dimension && "The dimension of the matrix was exceeded");

   // possibly add regularization
   if (this->use_regularization) {
      insert(ElementType(0), column_index, column_index);
   }
   this->current_column++;

   // start the next column at the current start
   if (column_index < this->dimension - 1) {
      this->column_starts[column_index + 2] = this->column_starts[column_index + 1];
   }
}

template <typename IndexType, typename ElementType>
ElementType CSCSymmetricMatrix<IndexType, ElementType>::smallest_diagonal_entry() const {
   ElementType smallest_entry = INF<ElementType>;
   for (size_t row_index: Range(this->dimension)) {
      smallest_entry = std::min(smallest_entry, this->diagonal_entries[row_index]);
   }
   return smallest_entry;
}

template <typename IndexType, typename ElementType>
void CSCSymmetricMatrix<IndexType, ElementType>::set_regularization(const std::function<ElementType(size_t /*index*/)>& regularization_function) {
   assert(this->use_regularization && "You are trying to regularize a matrix where regularization was not preallocated.");

   for (size_t row_index: Range(this->dimension)) {
      // the regularization term is located at the end of the column, that is right before the start of the next column
      const size_t k = this->column_starts[row_index + 1] - 1;
      const ElementType element = regularization_function(row_index);
      this->entries[k] = element;
      // update diagonal
      this->diagonal_entries[row_index] += element;
   }
}

template <typename ElementType>
std::tuple<size_t, size_t, ElementType> CSCSymmetricMatrix<ElementType>::dereference_iterator(size_t column_index, size_t nonzero_index) const {
   return {this->row_indices[nonzero_index], column_index, this->entries[nonzero_index]};
}

template <typename ElementType>
void CSCSymmetricMatrix<ElementType>::increment_iterator(size_t& column_index, size_t& nonzero_index) const {
   if (this->column_starts[column_index] <= nonzero_index && nonzero_index + 1 < this->column_starts[column_index + 1]) {
      // stay in the column
      nonzero_index++;
   }
   else {
      // move on to the next non-empty column
      do {
         column_index++;
      } while (column_index < this->dimension && this->column_starts[column_index] == this->column_starts[column_index + 1]);
      nonzero_index = this->column_starts[column_index];
   }
}

template <typename ElementType>
void CSCSymmetricMatrix<ElementType>::print(std::ostream& stream) const {
   stream << "W = "; print_vector(stream, this->entries, 0, this->number_nonzeros);
   stream << "with column start: "; print_vector(stream, this->column_starts, 0, this->dimension + 1);
   stream << "and row index: "; print_vector(stream, this->row_indices, 0, this->number_nonzeros);
}

template <typename IndexType, typename ElementType>
CSCSymmetricMatrix<IndexType, ElementType> CSCSymmetricMatrix<IndexType, ElementType>::identity(size_t dimension) {
   CSCSymmetricMatrix<IndexType, double> matrix(dimension, dimension, false);
   for (size_t row_index: Range(dimension)) {
      matrix.insert(1., row_index, row_index);
      matrix.finalize_column(row_index);
   }
   return matrix;
}

#endif // UNO_CSCSYMMETRICMATRIX_H
