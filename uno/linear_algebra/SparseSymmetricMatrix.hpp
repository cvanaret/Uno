// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SPARSESYMMETRICMATRIX_H
#define UNO_SPARSESYMMETRICMATRIX_H

#include <algorithm>
#include <memory>
#include <vector>
#include "SymmetricMatrix.hpp"

namespace uno {
   // forward declaration
   template <typename value_type>
   class Collection;

   // abstract class
   template <typename SparseStorage>
   class SparseSymmetricMatrix: public SymmetricMatrix<typename SparseStorage::index_type, typename SparseStorage::value_type> {
   public:
      using index_type = typename SparseStorage::index_type;
      using value_type = typename SparseStorage::value_type;

      SparseSymmetricMatrix(size_t dimension, size_t capacity, size_t regularization_size):
         SymmetricMatrix<index_type, value_type>(),
         sparse_storage(dimension, capacity, regularization_size) {
      }

      SparseSymmetricMatrix() = default;
      ~SparseSymmetricMatrix() override = default;
      SparseSymmetricMatrix& operator=(const SparseSymmetricMatrix& other) = default;
      SparseSymmetricMatrix& operator=(SparseSymmetricMatrix&& other) = default;

      void reset() override {
         this->sparse_storage.reset();
      }

      [[nodiscard]] size_t dimension() const override {
         return this->sparse_storage.dimension;
      }

      void set_dimension(size_t new_dimension) override {
         this->sparse_storage.set_dimension(new_dimension);
      }

      [[nodiscard]] size_t number_nonzeros() const override {
         return this->sparse_storage.number_nonzeros;
      }

      [[nodiscard]] size_t capacity() const override {
         return this->sparse_storage.capacity;
      }

      // build the matrix incrementally
      void insert(value_type term, index_type row_index, index_type column_index) override {
         this->sparse_storage.insert(term, row_index, column_index);
      }

      void finalize_column(index_type column_index) override {
         this->sparse_storage.finalize_column(column_index);
      }

      [[nodiscard]] value_type smallest_diagonal_entry(size_t max_dimension) const override {
         // diagonal entries might be at several locations and must be accumulated
         // TODO preallocate this vector somewhere
         std::vector<value_type> diagonal_entries(max_dimension, value_type(0));

         for (const auto [row_index, column_index, element]: *this) {
            if (row_index == column_index && row_index < max_dimension) {
               diagonal_entries[row_index] += element;
            }
         }
         // look at the first max_dimension diagonal entries
         return *std::min_element(diagonal_entries.begin(), diagonal_entries.end());
      }

      void set_regularization(const Collection<size_t>& indices, size_t offset, double factor) override {
         this->sparse_storage.set_regularization(indices, offset, factor);
      }

      [[nodiscard]] const value_type* data_pointer() const noexcept override { return this->sparse_storage.data_pointer(); }
      [[nodiscard]] value_type* data_pointer() noexcept override { return this->sparse_storage.data_pointer(); }

      void print(std::ostream& stream) const override { this->sparse_storage.print(stream); }

   protected:
      SparseStorage sparse_storage{};

      // virtual iterator functions
      [[nodiscard]] std::tuple<index_type, index_type, value_type> dereference_iterator(size_t column_index,
            size_t nonzero_index) const override {
         return this->sparse_storage.dereference_iterator(column_index, nonzero_index);
      }

      void increment_iterator(size_t& column_index, size_t& nonzero_index) const override {
         this->sparse_storage.increment_iterator(column_index, nonzero_index);
      }
   };
} // namespace

#endif // UNO_SPARSESYMMETRICMATRIX_H