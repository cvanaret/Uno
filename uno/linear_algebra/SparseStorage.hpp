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
            return a.nonzero_index != b.nonzero_index;
         }

      protected:
         const SparseStorage<IndexType, ElementType>& sparse_storage;
         size_t column_index;
         size_t nonzero_index;
      };

      using value_type = ElementType;

      size_t dimension;
      size_t number_nonzeros{0};
      size_t capacity;

      SparseStorage(size_t dimension, size_t capacity, bool use_regularization);
      virtual ~SparseStorage() = default;

      virtual void reset() = 0;
      void set_dimension(size_t new_dimension);

      // build the matrix incrementally
      virtual void insert(ElementType term, IndexType row_index, IndexType column_index) = 0;
      // this method will be used by the CSCSparseStorage subclass
      virtual void finalize_column(IndexType column_index) = 0;
      virtual void set_regularization(const std::function<ElementType(size_t /*index*/)>& regularization_function) = 0;
      virtual const ElementType* data_pointer() const noexcept = 0;
      virtual ElementType* data_pointer() noexcept = 0;

      [[nodiscard]] iterator begin() const {
         return iterator(*this, 0, 0);
      }
      [[nodiscard]] iterator end() const {
         return iterator(*this, this->dimension, this->number_nonzeros);
      }

      virtual void print(std::ostream& stream) const = 0;
      template <typename Index, typename Element>
      friend std::ostream& operator<<(std::ostream& stream, const SparseStorage<Index, Element>& matrix);

   protected:
      // regularization
      const bool use_regularization;

      // virtual iterator functions
      [[nodiscard]] virtual std::tuple<IndexType, IndexType, ElementType> dereference_iterator(size_t column_index, size_t nonzero_index) const = 0;
      virtual void increment_iterator(size_t& column_index, size_t& nonzero_index) const = 0;
   };

   // implementation

   template <typename IndexType, typename ElementType>
   SparseStorage<IndexType, ElementType>::SparseStorage(size_t dimension, size_t capacity, bool use_regularization) :
         dimension(dimension),
         // if regularization is used, allocate the necessary space
         capacity(capacity + (use_regularization ? dimension : 0)),
         use_regularization(use_regularization) {
   }

   template <typename IndexType, typename ElementType>
   void SparseStorage<IndexType, ElementType>::set_dimension(size_t new_dimension) {
      this->dimension = new_dimension;
   }
   
   template <typename Index, typename Element>
   std::ostream& operator<<(std::ostream& stream, const SparseStorage<Index, Element>& matrix) {
      stream << "Dimension: " << matrix.dimension << ", number of nonzeros: " << matrix.number_nonzeros << '\n';
      matrix.print(stream);
      return stream;
   }
} // namespace

#endif // UNO_SPARSESTORAGE_H
