// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SPARSESTORAGE_H
#define UNO_SPARSESTORAGE_H

#include <ostream>

namespace uno {
   // forward declaration
   template <typename ElementType>
   class Collection;

   // abstract class
   template <typename IndexType, typename ElementType>
   class SparseStorage {
   public:
      using index_type = IndexType;
      using value_type = ElementType;

      size_t dimension{0};
      size_t number_nonzeros{0};
      size_t capacity{0};

      SparseStorage(size_t dimension, size_t capacity, size_t regularization_size);
      SparseStorage() = default;
      SparseStorage& operator=(const SparseStorage& other) = default;
      SparseStorage& operator=(SparseStorage&& other) = default;
      virtual ~SparseStorage() = default;

      virtual void reset() = 0;
      virtual void set_dimension(size_t new_dimension) = 0;

      // build the matrix incrementally
      virtual void insert(IndexType row_index, IndexType column_index, ElementType term) = 0;
      // this method will be used by the CSCSparseStorage subclass
      virtual void finalize_column(IndexType column_index) = 0;
      virtual void set_regularization(const Collection<size_t>& indices, size_t offset, double factor) = 0;
      virtual const ElementType* data_pointer() const noexcept = 0;
      virtual ElementType* data_pointer() noexcept = 0;
      virtual void print(std::ostream& stream) const = 0;

      // virtual iterator functions
      [[nodiscard]] virtual std::tuple<IndexType, IndexType, ElementType> dereference_iterator(size_t column_index, size_t nonzero_index) const = 0;
      virtual void increment_iterator(size_t& column_index, size_t& nonzero_index) const = 0;

      template <typename Index, typename Element>
      friend std::ostream& operator<<(std::ostream& stream, const SparseStorage<Index, Element>& sparse_storage);

   protected:
      size_t regularization_size{0};
   };

   // implementation

   template <typename IndexType, typename ElementType>
   SparseStorage<IndexType, ElementType>::SparseStorage(size_t dimension, size_t capacity, size_t regularization_size) :
         dimension(dimension),
         // if regularization is used, allocate the necessary space
         capacity(capacity + regularization_size),
         regularization_size(regularization_size) {
   }
} // namespace

#endif // UNO_SPARSESTORAGE_H