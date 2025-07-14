// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_COOSPARSESTORAGE_H
#define UNO_COOSPARSESTORAGE_H

#include <cassert>
#include <vector>
#include "SparseStorage.hpp"
#include "symbolic/Range.hpp"
#include "symbolic/VectorView.hpp"

namespace uno {
   /*
    * Coordinate list
    * https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)
    */
   template <typename IndexType, typename ElementType>
   class COOFormat : public SparseStorage<IndexType, ElementType> {
   public:
      COOFormat(size_t dimension, size_t capacity, size_t regularization_size);
      COOFormat() = default;
      ~COOFormat() override = default;
      COOFormat& operator=(const COOFormat& other) = default;
      COOFormat& operator=(COOFormat&& other) = default;

      void reset() override;
      void set_dimension(size_t new_dimension) override;

      void insert(IndexType row_index, IndexType column_index, ElementType term) override;
      void finalize_column(IndexType /*column_index*/) override { /* do nothing */ }
      void set_regularization(const Collection<size_t>& indices, size_t offset, double factor) override;
      const ElementType* data_pointer() const noexcept override { return this->entries.data(); }
      ElementType* data_pointer() noexcept override { return this->entries.data(); }

      void print(std::ostream& stream) const override;

      // iterator functions
      [[nodiscard]] std::tuple<IndexType, IndexType, ElementType> dereference_iterator(size_t column_index, size_t nonzero_index) const override;
      void increment_iterator(size_t& column_index, size_t& nonzero_index) const override;

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
      std::vector<ElementType> entries{};
      std::vector<IndexType> row_indices{};
      std::vector<IndexType> column_indices{};

      void initialize_regularization();
   };

   // implementation

   template <typename IndexType, typename ElementType>
   COOFormat<IndexType, ElementType>::COOFormat(size_t dimension, size_t capacity, size_t regularization_size):
         SparseStorage<IndexType, ElementType>(dimension, capacity, regularization_size) {
      this->entries.reserve(this->capacity);
      this->row_indices.reserve(this->capacity);
      this->column_indices.reserve(this->capacity);

      // initialize regularization terms
      if (0 < this->regularization_size) {
         this->initialize_regularization();
      }
   }

   template <typename IndexType, typename ElementType>
   void COOFormat<IndexType, ElementType>::reset() {
      // empty the matrix
      this->number_nonzeros = 0;
      this->entries.clear();
      this->row_indices.clear();
      this->column_indices.clear();

      // initialize regularization terms
      if (0 < this->regularization_size) {
         this->initialize_regularization();
      }
   }

   template <typename IndexType, typename ElementType>
   void COOFormat<IndexType, ElementType>::set_dimension(size_t new_dimension) {
      this->dimension = new_dimension;
      this->entries.reserve(this->dimension);
      this->row_indices.reserve(this->dimension);
      this->column_indices.reserve(this->dimension);
   }

   template <typename IndexType, typename ElementType>
   void COOFormat<IndexType, ElementType>::insert(IndexType row_index, IndexType column_index, ElementType term) {
      assert(this->number_nonzeros <= this->row_indices.size() && "The COO matrix doesn't have a sufficient capacity");

      this->entries.emplace_back(term);
      this->row_indices.emplace_back(row_index);
      this->column_indices.emplace_back(column_index);
      this->number_nonzeros++;
   }

   template <typename IndexType, typename ElementType>
   void COOFormat<IndexType, ElementType>::set_regularization(const Collection<size_t>& indices, size_t offset, double factor) {
      assert(0 < this->regularization_size && "You are trying to regularize a matrix where regularization was not preallocated.");

      // the regularization terms (that lie at the start of the entries vector) can be directly modified
      for (size_t index: indices) {
         const ElementType element = factor;
         this->entries[index + offset] = element;
      }
   }

   template <typename IndexType, typename ElementType>
   void COOFormat<IndexType, ElementType>::print(std::ostream& stream) const {
      stream << "values = "; print_vector(stream, view(this->entries, 0, this->number_nonzeros));
      stream << "with column indices: "; print_vector(stream, view(this->row_indices, 0, this->number_nonzeros));
      stream << "and row indices: "; print_vector(stream, view(this->column_indices, 0, this->number_nonzeros));
   }

   template <typename IndexType, typename ElementType>
   void COOFormat<IndexType, ElementType>::initialize_regularization() {
      // introduce elements at the start of the entries
      for (size_t row_index: Range(this->regularization_size)) {
         this->insert(IndexType(row_index), IndexType(row_index), ElementType(0));
      }
   }

   template <typename IndexType, typename ElementType>
   std::tuple<IndexType, IndexType, ElementType> COOFormat<IndexType, ElementType>::dereference_iterator(size_t /*column_index*/,
         size_t nonzero_index) const {
      return {this->row_indices[nonzero_index], this->column_indices[nonzero_index], this->entries[nonzero_index]};
   }

   template <typename IndexType, typename ElementType>
   void COOFormat<IndexType, ElementType>::increment_iterator(size_t& column_index, size_t& nonzero_index) const {
      nonzero_index++;
      // if end reached
      if (nonzero_index == this->number_nonzeros) {
         column_index = this->dimension;
      }
   }
} // namespace

#endif // UNO_COOSPARSESTORAGE_H