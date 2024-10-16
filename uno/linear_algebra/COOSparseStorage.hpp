// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_COOSPARSESTORAGE_H
#define UNO_COOSPARSESTORAGE_H

#include <cassert>
#include <vector>
#include "SparseStorage.hpp"
#include "symbolic/Range.hpp"

namespace uno {
   /*
    * Coordinate list
    * https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)
    */
   template <typename IndexType, typename ElementType>
   class COOSparseStorage : public SparseStorage<IndexType, ElementType> {
   public:
      COOSparseStorage(size_t dimension, size_t capacity, bool use_regularization);

      void reset() override;
      void insert(ElementType term, IndexType row_index, IndexType column_index) override;
      void finalize_column(IndexType /*column_index*/) override { /* do nothing */ }
      void set_regularization(const std::function<ElementType(size_t index)>& regularization_function) override;
      const ElementType* data_pointer() const noexcept override { return this->entries.data(); }
      ElementType* data_pointer() noexcept override { return this->entries.data(); }

      void print(std::ostream& stream) const override;

      const IndexType* row_indices_pointer() const {
         return this->row_indices.data();
      }
      IndexType* row_indices_pointer() {
         return this->row_indices.data();
      }
      const IndexType* column_indices_pointer() const {
         return this->column_indices.data();
      }
      IndexType* column_indices_pointer() {
         return this->column_indices.data();
      }

   protected:
      std::vector<ElementType> entries;
      std::vector<IndexType> row_indices;
      std::vector<IndexType> column_indices;

      void initialize_regularization();

      // iterator functions
      [[nodiscard]] std::tuple<IndexType, IndexType, ElementType> dereference_iterator(size_t column_index, size_t nonzero_index) const override;
      void increment_iterator(size_t& column_index, size_t& nonzero_index) const override;
   };

   // implementation

   template <typename IndexType, typename ElementType>
   COOSparseStorage<IndexType, ElementType>::COOSparseStorage(size_t dimension, size_t capacity, bool use_regularization):
         SparseStorage<IndexType, ElementType>(dimension, capacity, use_regularization) {
      this->entries.reserve(this->capacity);
      this->row_indices.reserve(this->capacity);
      this->column_indices.reserve(this->capacity);

      // initialize regularization terms
      if (this->use_regularization) {
         this->initialize_regularization();
      }
   }

   template <typename IndexType, typename ElementType>
   void COOSparseStorage<IndexType, ElementType>::reset() {
      // empty the matrix
      this->number_nonzeros = 0;
      this->entries.clear();
      this->row_indices.clear();
      this->column_indices.clear();

      // initialize regularization terms
      if (this->use_regularization) {
         this->initialize_regularization();
      }
   }

   template <typename IndexType, typename ElementType>
   void COOSparseStorage<IndexType, ElementType>::insert(ElementType term, IndexType row_index, IndexType column_index) {
      assert(this->number_nonzeros <= this->row_indices.size() && "The COO matrix doesn't have a sufficient capacity");

      this->entries.emplace_back(term);
      this->row_indices.emplace_back(row_index);
      this->column_indices.emplace_back(column_index);
      this->number_nonzeros++;
   }

   template <typename IndexType, typename ElementType>
   void COOSparseStorage<IndexType, ElementType>::set_regularization(const std::function<ElementType(size_t /*index*/)>& regularization_function) {
      assert(this->use_regularization && "You are trying to regularize a matrix where regularization was not preallocated.");

      // the regularization terms (that lie at the start of the entries vector) can be directly modified
      for (size_t row_index: Range(this->dimension)) {
         const ElementType element = regularization_function(row_index);
         this->entries[row_index] = element;
      }
   }

   template <typename IndexType, typename ElementType>
   void COOSparseStorage<IndexType, ElementType>::print(std::ostream& stream) const {
      for (const auto [row_index, column_index, element]: *this) {
         stream << "m(" << row_index << ", " << column_index << ") = " << element << '\n';
      }
   }

   template <typename IndexType, typename ElementType>
   void COOSparseStorage<IndexType, ElementType>::initialize_regularization() {
      // introduce elements at the start of the entries
      for (size_t row_index: Range(this->dimension)) {
         this->insert(ElementType(0), IndexType(row_index), IndexType(row_index));
      }
   }

   template <typename IndexType, typename ElementType>
   std::tuple<IndexType, IndexType, ElementType> COOSparseStorage<IndexType, ElementType>::dereference_iterator(size_t /*column_index*/,
         size_t nonzero_index) const {
      return {this->row_indices[nonzero_index], this->column_indices[nonzero_index], this->entries[nonzero_index]};
   }

   template <typename IndexType, typename ElementType>
   void COOSparseStorage<IndexType, ElementType>::increment_iterator(size_t& column_index, size_t& nonzero_index) const {
      nonzero_index++;
      // if end reached
      if (nonzero_index == this->number_nonzeros) {
         column_index = this->dimension;
      }
   }
} // namespace

#endif // UNO_COOSPARSESTORAGE_H
