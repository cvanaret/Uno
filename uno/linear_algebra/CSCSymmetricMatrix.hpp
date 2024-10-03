// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_CSCSYMMETRICMATRIX_H
#define UNO_CSCSYMMETRICMATRIX_H

#include <ostream>
#include <cassert>
#include "SymmetricMatrix.hpp"
#include "tools/Infinity.hpp"
#include "symbolic/VectorView.hpp"

namespace uno {
   /*
 * Compressed Sparse Column
 * https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
 */

   template <typename IndexType, typename ElementType>
   class CSCSymmetricMatrix : public SymmetricMatrix<IndexType, ElementType> {
   public:
      CSCSymmetricMatrix(size_t dimension, size_t original_capacity, bool use_regularization);

      void reset() override;
      void insert(ElementType term, IndexType row_index, IndexType column_index) override;
      void finalize_column(IndexType column_index) override;
      void set_regularization(const std::function<ElementType(IndexType /*index*/)>& regularization_function) override;

      void print(std::ostream& stream) const override;

      static CSCSymmetricMatrix<IndexType, ElementType> identity(size_t dimension);

   protected:
      // entries and row_indices have nnz elements
      // column_starts has dimension+1 elements
      Vector<IndexType> column_starts{};
      std::vector<IndexType> row_indices{};
      IndexType current_column{0};

      // iterator functions
      [[nodiscard]] std::tuple<IndexType, IndexType, ElementType> dereference_iterator(IndexType column_index, size_t nonzero_index) const override;
      void increment_iterator(IndexType& column_index, IndexType& nonzero_index) const override;
   };

   template <typename IndexType, typename ElementType>
   CSCSymmetricMatrix<IndexType, ElementType>::CSCSymmetricMatrix(size_t dimension, size_t original_capacity, bool use_regularization):
         SymmetricMatrix<IndexType, ElementType>(dimension, original_capacity, use_regularization),
         column_starts(dimension + 1) {
      this->entries.reserve(this->capacity);
      this->row_indices.reserve(this->capacity);
   }

   template <typename IndexType, typename ElementType>
   void CSCSymmetricMatrix<IndexType, ElementType>::reset() {
      // empty the matrix
      SymmetricMatrix<IndexType, ElementType>::reset();
      this->row_indices.clear();
      this->column_starts.fill(0);
      this->current_column = 0;
   }

   template <typename IndexType, typename ElementType>
   void CSCSymmetricMatrix<IndexType, ElementType>::insert(ElementType term, IndexType row_index, IndexType column_index) {
      assert(column_index == this->current_column && "The previous columns should be finalized");

      this->entries.emplace_back(term);
      this->row_indices.emplace_back(row_index);
      this->column_starts[column_index + 1]++;
      this->number_nonzeros++;
   }

   template <typename IndexType, typename ElementType>
   void CSCSymmetricMatrix<IndexType, ElementType>::finalize_column(IndexType column_index) {
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
   void CSCSymmetricMatrix<IndexType, ElementType>::set_regularization(const std::function<ElementType(IndexType /*index*/)>& regularization_function) {
      assert(this->use_regularization && "You are trying to regularize a matrix where regularization was not preallocated.");

      for (size_t row_index: Range(this->dimension)) {
         // the regularization term is located at the end of the column, that is right before the start of the next column
         const size_t k = static_cast<size_t>(this->column_starts[row_index + 1] - 1);
         const ElementType element = regularization_function(row_index);
         this->entries[k] = element;
      }
   }

   template <typename IndexType, typename ElementType>
   std::tuple<IndexType, IndexType, ElementType> CSCSymmetricMatrix<IndexType, ElementType>::dereference_iterator(IndexType column_index,
         size_t nonzero_index) const {
      return {this->row_indices[nonzero_index], column_index, this->entries[nonzero_index]};
   }

   template <typename IndexType, typename ElementType>
   void CSCSymmetricMatrix<IndexType, ElementType>::increment_iterator(IndexType& column_index, IndexType& nonzero_index) const {
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

   template <typename IndexType, typename ElementType>
   void CSCSymmetricMatrix<IndexType, ElementType>::print(std::ostream& stream) const {
      stream << "W = "; print_vector(stream, view(this->entries, 0, this->number_nonzeros));
      stream << "with column start: "; print_vector(stream, view(this->column_starts, 0, this->dimension + 1));
      stream << "and row index: "; print_vector(stream, view(this->row_indices, 0, this->number_nonzeros));
   }

   template <typename IndexType, typename ElementType>
   CSCSymmetricMatrix<IndexType, ElementType> CSCSymmetricMatrix<IndexType, ElementType>::identity(size_t dimension) {
      CSCSymmetricMatrix<IndexType, ElementType> matrix(dimension, dimension, false);
      for (size_t row_index: Range(dimension)) {
         matrix.insert(ElementType(1), IndexType(row_index), IndexType(row_index));
         matrix.finalize_column(row_index);
      }
      return matrix;
   }
} // namespace

#endif // UNO_CSCSYMMETRICMATRIX_H