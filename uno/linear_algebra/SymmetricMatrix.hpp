// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICMATRIX_H
#define UNO_SYMMETRICMATRIX_H

#include <cassert>
#include <memory>
#include "SparseStorageFactory.hpp"

namespace uno {
   // forward declaration
   template <typename ElementType>
   class Collection;

   // abstract class
   template <typename IndexType, typename ElementType>
   class SymmetricMatrix {
   public:
      using value_type = ElementType;

      class iterator {
      public:
         iterator(const SymmetricMatrix& matrix, size_t column_index, size_t nonzero_index) :
               matrix(matrix), column_index(column_index), nonzero_index(nonzero_index) {
         }

         [[nodiscard]] std::tuple<IndexType, IndexType, ElementType> operator*() const {
            return this->matrix.dereference_iterator(this->column_index, this->nonzero_index);
         }

         iterator& operator++() {
            this->matrix.increment_iterator(this->column_index, this->nonzero_index);
            return *this;
         }

         friend bool operator!=(const iterator& a, const iterator& b) {
            return a.nonzero_index != b.nonzero_index;
         }

      protected:
         const SymmetricMatrix& matrix;
         size_t column_index;
         size_t nonzero_index;
      };

      SymmetricMatrix() = default;

      virtual ~SymmetricMatrix() = default;
      SymmetricMatrix& operator=(SymmetricMatrix&& other) = default;

      virtual void reset() = 0;
      [[nodiscard]] virtual size_t dimension() const = 0;
      virtual void set_dimension(size_t new_dimension) = 0;
      [[nodiscard]] virtual size_t number_nonzeros() const = 0;
      [[nodiscard]] virtual size_t capacity() const = 0;

      // build the matrix incrementally
      virtual void insert(ElementType term, IndexType row_index, IndexType column_index) = 0;
      virtual void finalize_column(IndexType column_index) = 0;
      
      [[nodiscard]] virtual ElementType smallest_diagonal_entry(size_t max_dimension) const = 0;
      
      virtual void set_regularization(const Collection<size_t>& indices, size_t offset, double factor) = 0;

      iterator begin() const {
         return iterator(*this, 0, 0);
      }

      iterator end() const {
         return iterator(*this, this->dimension(), this->number_nonzeros());
      }

      [[nodiscard]] virtual const ElementType* data_pointer() const noexcept = 0;
      [[nodiscard]] virtual ElementType* data_pointer() noexcept = 0;

      template <typename Vector1, typename Vector2>
      void product(const Vector1& vector, Vector2& result) const;
      template <typename Vector1, typename Vector2>
      ElementType quadratic_product(const Vector1& x, const Vector2& y) const;

      virtual void print(std::ostream& stream) const = 0;
      template <typename Index, typename Element>
      friend std::ostream& operator<<(std::ostream& stream, const SymmetricMatrix<Index, Element>& matrix);

   protected:
      // virtual iterator functions
      [[nodiscard]] virtual std::tuple<IndexType, IndexType, ElementType> dereference_iterator(size_t column_index,
         size_t nonzero_index) const = 0;
      virtual void increment_iterator(size_t& column_index, size_t& nonzero_index) const = 0;
   };

   // implementation

   template <typename IndexType, typename ElementType>
   template <typename Vector1, typename Vector2>
   inline void SymmetricMatrix<IndexType, ElementType>::product(const Vector1& vector, Vector2& result) const {
      for (const auto [row_index, column_index, element]: *this) { // get the (row_index, column_index) element
         result[row_index] += element * vector[column_index];
         if (row_index != column_index) {
            // consider the (column_index, row_index) as well (by symmetry)
            result[column_index] += element * vector[row_index];
         }
      }
   }

   template <typename IndexType, typename ElementType>
   template <typename Vector1, typename Vector2>
   inline ElementType SymmetricMatrix<IndexType, ElementType>::quadratic_product(const Vector1& x, const Vector2& y) const {
      static_assert(std::is_same_v<typename Vector1::value_type, ElementType>);
      static_assert(std::is_same_v<typename Vector2::value_type, ElementType>);
      assert(x.size() == y.size() && "SymmetricMatrix::quadratic_product: the two vectors x and y do not have the same size");

      ElementType result = ElementType(0);
      for (const auto [row_index, column_index, element]: *this) {
         if (row_index == column_index) {
            // diagonal term
            result += element * x[row_index] * y[row_index];
         }
         else {
            // off-diagonal term
            result += element * (x[row_index] * y[column_index] + x[column_index] * y[row_index]);
         }
      }
      return result;
   }

   template <typename IndexType, typename ElementType>
   inline std::ostream& operator<<(std::ostream& stream, const SymmetricMatrix<IndexType, ElementType>& matrix) {
      stream << "Dimension: " << matrix.dimension() << ", number of nonzeros: " << matrix.number_nonzeros() << '\n';
      matrix.print(stream);
      return stream;
   }
} // namespace

#endif // UNO_SYMMETRICMATRIX_H