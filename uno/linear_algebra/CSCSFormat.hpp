// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_CSCSPARSESTORAGE_H
#define UNO_CSCSPARSESTORAGE_H

#include <cassert>
#include <vector>
#include "SparseStorage.hpp"
#include "linear_algebra/Vector.hpp"
#include "symbolic/VectorView.hpp"

namespace uno {
   /*
    * Compressed Sparse Column
    * https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
    */
   template <typename IndexType, typename ElementType>
   class CSCSFormat : public SparseStorage<IndexType, ElementType> {
   public:
      CSCSFormat(size_t dimension, size_t capacity, size_t regularization_size);
      CSCSFormat() = default;
      ~CSCSFormat() override = default;
      CSCSFormat& operator=(const CSCSFormat& other) = default;
      CSCSFormat& operator=(CSCSFormat&& other) = default;

      void reset() override;
      void set_dimension(size_t new_dimension) override;

      void insert(ElementType term, IndexType row_index, IndexType column_index) override;
      void finalize_column(IndexType column_index) override;
      void set_regularization(const Collection<size_t>& indices, size_t offset, double factor) override;
      const ElementType* data_pointer() const noexcept override { return this->entries.data(); }
      ElementType* data_pointer() noexcept override { return this->entries.data(); }

      void print(std::ostream& stream) const override;

      // iterator functions
      [[nodiscard]] std::tuple<IndexType, IndexType, ElementType> dereference_iterator(IndexType column_index, size_t nonzero_index) const override;
      void increment_iterator(IndexType& column_index, IndexType& nonzero_index) const override;

   protected:
      std::vector<ElementType> entries{};
      // entries and row_indices have nnz elements
      // column_starts has dimension+1 elements
      Vector<IndexType> column_starts{};
      std::vector<IndexType> row_indices{};
      IndexType current_column{0};
   };

   template <typename IndexType, typename ElementType>
   CSCSFormat<IndexType, ElementType>::CSCSFormat(size_t dimension, size_t capacity, size_t regularization_size):
         SparseStorage<IndexType, ElementType>(dimension, capacity, regularization_size),
         column_starts(dimension + 1) {
      this->entries.reserve(this->capacity);
      this->row_indices.reserve(this->capacity);
   }

   template <typename IndexType, typename ElementType>
   void CSCSFormat<IndexType, ElementType>::reset() {
      // empty the matrix
      this->number_nonzeros = 0;
      this->entries.clear();
      this->row_indices.clear();
      this->column_starts.fill(0);
      this->current_column = 0;
   }

   template <typename IndexType, typename ElementType>
   void CSCSFormat<IndexType, ElementType>::set_dimension(size_t new_dimension) {
      this->column_starts.resize(new_dimension + 1);
      // if we enlarge the storage, copy the column starts
      if (this->dimension < new_dimension) {
         for (size_t index: Range(this->dimension + 1, new_dimension + 1)) {
            this->column_starts[index] = this->column_starts[this->dimension];
         }
      }
      this->dimension = new_dimension;
   }

   template <typename IndexType, typename ElementType>
   void CSCSFormat<IndexType, ElementType>::insert(ElementType term, IndexType row_index, IndexType column_index) {
      assert(column_index == this->current_column && "The previous columns should be finalized");

      this->entries.emplace_back(term);
      this->row_indices.emplace_back(row_index);
      this->column_starts[column_index + 1]++;
      this->number_nonzeros++;
   }

   template <typename IndexType, typename ElementType>
   void CSCSFormat<IndexType, ElementType>::finalize_column(IndexType column_index) {
      assert(column_index == this->current_column && "You are not finalizing the current column");
      assert(column_index < this->dimension && "The dimension of the matrix was exceeded");

      // possibly add regularization
      if (0 < this->regularization_size) {
         this->insert(ElementType(0), column_index, column_index);
      }
      this->current_column++;

      // start the next column at the current start
      if (column_index < this->dimension - 1) {
         this->column_starts[column_index + 2] = this->column_starts[column_index + 1];
      }
   }

   template <typename IndexType, typename ElementType>
   void CSCSFormat<IndexType, ElementType>::set_regularization(const Collection<size_t>& indices, size_t offset, double factor) {
      assert(0 < this->regularization_size && "You are trying to regularize a matrix where regularization was not preallocated.");

      for (size_t row_index: indices) {
         // the regularization term is located at the end of the column, that is right before the start of the next column
         const size_t position = static_cast<size_t>(this->column_starts[row_index + offset + 1] - 1);
         const ElementType element = factor;
         this->entries[position] = element;
      }
   }

   template <typename IndexType, typename ElementType>
   std::tuple<IndexType, IndexType, ElementType> CSCSFormat<IndexType, ElementType>::dereference_iterator(IndexType column_index,
         size_t nonzero_index) const {
      return {this->row_indices[nonzero_index], column_index, this->entries[nonzero_index]};
   }

   template <typename IndexType, typename ElementType>
   void CSCSFormat<IndexType, ElementType>::increment_iterator(IndexType& column_index, IndexType& nonzero_index) const {
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
   void CSCSFormat<IndexType, ElementType>::print(std::ostream& stream) const {
      stream << "W = "; print_vector(stream, view(this->entries, 0, this->number_nonzeros));
      stream << "with column start: "; print_vector(stream, view(this->column_starts, 0, this->dimension + 1));
      stream << "and row index: "; print_vector(stream, view(this->row_indices, 0, this->number_nonzeros));
   }
} // namespace

#endif // UNO_CSCSPARSESTORAGE_H