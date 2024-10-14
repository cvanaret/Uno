// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FORTRANCOOSPARSESTORAGE_H
#define UNO_FORTRANCOOSPARSESTORAGE_H

#include <cassert>
#include "linear_algebra/SparseStorage.hpp"
#include "symbolic/Range.hpp"

namespace uno {
   /*
    * Coordinate list
    * https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)
    * indices are integers (int) and start at 1 (Fortran convention)
    */
   template <typename IndexType, typename ElementType>
   class FortranCOOSparseStorage : public SparseStorage<IndexType, ElementType> {
   public:
      FortranCOOSparseStorage(size_t number_rows, size_t number_columns, size_t capacity, bool use_regularization);

      void reset() override;
      void insert(ElementType term, IndexType row_index, IndexType column_index) override;
      void finalize_column(IndexType /*column_index*/) override { /* do nothing */ }
      void set_regularization(const std::function<ElementType(IndexType index)>& regularization_function) override;
      const ElementType* data_pointer() const noexcept override { return this->entries.data(); }
      ElementType* data_pointer() noexcept override { return this->entries.data(); }

      void print(std::ostream& stream) const override;

      [[nodiscard]] const int* row_indices_pointer() const {
         return this->row_indices.data();
      }
      int* row_indices_pointer() {
         return this->row_indices.data();
      }
      [[nodiscard]] const int* column_indices_pointer() const {
         return this->column_indices.data();
      }
      int* column_indices_pointer() {
         return this->column_indices.data();
      }

   protected:
      std::vector<ElementType> entries;
      std::vector<int> row_indices;
      std::vector<int> column_indices;
      const int fortran_shift{1};

      void initialize_regularization();

      // iterator functions
      [[nodiscard]] std::tuple<IndexType, IndexType, ElementType> dereference_iterator(size_t column_index, size_t nonzero_index) const override;
      void increment_iterator(size_t& column_index, size_t& nonzero_index) const override;
   };

   // implementation

   template <typename IndexType, typename ElementType>
   FortranCOOSparseStorage<IndexType, ElementType>::FortranCOOSparseStorage(size_t number_rows, size_t number_columns, size_t capacity,
         bool use_regularization): SparseStorage<IndexType, ElementType>(number_rows, number_columns, capacity, use_regularization) {
      assert((not use_regularization || number_rows == number_columns) && "FortranCOOSparseStorage: a non-square matrix cannot be regularized.");

      this->entries.reserve(this->capacity);
      this->row_indices.reserve(this->capacity);
      this->column_indices.reserve(this->capacity);

      // initialize regularization terms
      if (this->use_regularization) {
         this->initialize_regularization();
      }
   }

   template <typename IndexType, typename ElementType>
   void FortranCOOSparseStorage<IndexType, ElementType>::reset() {
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
   void FortranCOOSparseStorage<IndexType, ElementType>::insert(ElementType term, IndexType row_index, IndexType column_index) {
      assert(this->number_nonzeros <= this->row_indices.size() && "The COO matrix doesn't have a sufficient capacity");

      this->entries.emplace_back(term);
      this->row_indices.emplace_back(static_cast<int>(row_index) + this->fortran_shift);
      this->column_indices.emplace_back(static_cast<int>(column_index) + this->fortran_shift);
      this->number_nonzeros++;
   }

   template <typename IndexType, typename ElementType>
   void FortranCOOSparseStorage<IndexType, ElementType>::set_regularization(const std::function<ElementType(IndexType /*index*/)>& regularization_function) {
      assert(this->use_regularization && "You are trying to regularize a matrix where regularization was not preallocated.");

      // the regularization terms (that lie at the start of the entries vector) can be directly modified
      for (size_t row_index: Range(std::min(this->number_rows, this->number_columns))) {
         const ElementType element = regularization_function(static_cast<IndexType>(row_index));
         this->entries[row_index] = element;
      }
   }

   template <typename IndexType, typename ElementType>
   void FortranCOOSparseStorage<IndexType, ElementType>::print(std::ostream& stream) const {
      for (size_t index = 0; index < this->number_nonzeros; index++) {
         stream << "m(" << this->row_indices[index] << ", " << this->column_indices[index] << ") = " << this->entries[index] << '\n';
      }
   }

   template <typename IndexType, typename ElementType>
   void FortranCOOSparseStorage<IndexType, ElementType>::initialize_regularization() {
      // introduce elements at the start of the entries
      for (size_t row_index: Range(std::min(this->number_rows, this->number_columns))) {
         this->insert(ElementType(0), IndexType(row_index), IndexType(row_index));
      }
   }

   template <typename IndexType, typename ElementType>
   std::tuple<IndexType, IndexType, ElementType> FortranCOOSparseStorage<IndexType, ElementType>::dereference_iterator(size_t /*column_index*/,
         size_t nonzero_index) const {
      return {static_cast<IndexType>(this->row_indices[nonzero_index] - this->fortran_shift),
            static_cast<IndexType>(this->column_indices[nonzero_index] - this->fortran_shift),
            this->entries[nonzero_index]};
   }

   template <typename IndexType, typename ElementType>
   void FortranCOOSparseStorage<IndexType, ElementType>::increment_iterator(size_t& column_index, size_t& nonzero_index) const {
      nonzero_index++;
      // if end reached
      if (nonzero_index == this->number_nonzeros) {
         column_index = this->number_columns;
      }
   }
} // namespace

#endif // UNO_FORTRANCOOSPARSESTORAGE_H
