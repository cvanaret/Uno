// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SPARSESTORAGE_H
#define UNO_SPARSESTORAGE_H

#include <ostream>
#include <functional>

namespace uno {
   // abstract class
   template <typename IndexType, typename ElementType>
   class SparseStorage {
   public:
      class iterator {
      public:
         iterator(const SparseStorage<IndexType, ElementType>& sparse_storage, size_t column_index, size_t nonzero_index) :
               sparse_storage(sparse_storage), column_index(column_index), nonzero_index(nonzero_index) {
         }

         [[nodiscard]] std::tuple<IndexType, IndexType, ElementType> operator*() const {
            return this->sparse_storage.dereference_iterator(this->column_index, this->nonzero_index);
         }

         iterator& operator++() {
            this->sparse_storage.increment_iterator(this->column_index, this->nonzero_index);
            return *this;
         }

         friend bool operator!=(const iterator& a, const iterator& b) {
            return &a.sparse_storage != &b.sparse_storage || a.column_index != b.column_index || a.nonzero_index != b.nonzero_index;
         }

      protected:
         const SparseStorage<IndexType, ElementType>& sparse_storage;
         size_t column_index;
         size_t nonzero_index;
      };

      using value_type = ElementType;

      SparseStorage(size_t number_rows, size_t number_columns, size_t capacity, bool use_regularization);
      virtual ~SparseStorage() = default;

      virtual void reset() = 0;
      void set_number_rows(size_t new_number_rows) { this->number_rows = new_number_rows; }
      void set_number_columns(size_t new_number_columns) { this->number_columns = new_number_columns; }
      [[nodiscard]] size_t get_number_rows() const { return this->number_rows; }
      [[nodiscard]] size_t get_number_columns() const { return this->number_columns; }
      [[nodiscard]] size_t get_number_nonzeros() const { return this->number_nonzeros; }
      [[nodiscard]] size_t get_capacity() const { return this->capacity; }

      // build the matrix incrementally
      virtual void insert(ElementType term, IndexType row_index, IndexType column_index) = 0;
      // this method will be used by the CSCSparseStorage subclass
      virtual void finalize_column(IndexType column_index) = 0;
      virtual void set_regularization(const std::function<ElementType(IndexType /*index*/)>& regularization_function) = 0;
      virtual const ElementType* data_pointer() const noexcept = 0;
      virtual ElementType* data_pointer() noexcept = 0;

      [[nodiscard]] iterator begin() const {
         return iterator(*this, 0, 0);
      }
      [[nodiscard]] iterator end() const {
         return iterator(*this, this->number_columns, this->number_nonzeros);
      }

      virtual void print(std::ostream& stream) const = 0;
      template <typename Index, typename Element>
      friend std::ostream& operator<<(std::ostream& stream, const SparseStorage<Index, Element>& sparse_storage);

   protected:
      size_t number_rows, number_columns;
      size_t number_nonzeros{0};
      size_t capacity;
      const bool use_regularization;

      // virtual iterator functions
      [[nodiscard]] virtual std::tuple<IndexType, IndexType, ElementType> dereference_iterator(size_t column_index, size_t nonzero_index) const = 0;
      virtual void increment_iterator(size_t& column_index, size_t& nonzero_index) const = 0;
   };

   // implementation

   template <typename IndexType, typename ElementType>
   SparseStorage<IndexType, ElementType>::SparseStorage(size_t number_rows, size_t number_columns, size_t capacity, bool use_regularization) :
         number_rows(number_rows),
         number_columns(number_columns),
         // if regularization is used, allocate the necessary space
         capacity(capacity + (use_regularization ? std::min(number_rows, number_columns) : 0)),
         use_regularization(use_regularization) {
   }

   template <typename Index, typename Element>
   std::ostream& operator<<(std::ostream& stream, const SparseStorage<Index, Element>& sparse_storage) {
      stream << "Dimensions: (" << sparse_storage.number_rows << ", " << sparse_storage.number_columns << "), number of nonzeros: " << sparse_storage.number_nonzeros << '\n';
      sparse_storage.print(stream);
      return stream;
   }
} // namespace

#endif // UNO_SPARSESTORAGE_H
