// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICMATRIX_H
#define UNO_SYMMETRICMATRIX_H

#include <vector>
#include <functional>
#include <cassert>
#include "Matrix.hpp"
#include "Vector.hpp"
#include "SparseVector.hpp"

// abstract class
template <typename IndexType, typename ElementType>
class SymmetricMatrix: public Matrix<IndexType, ElementType> {
public:
   class iterator {
   public:
      iterator(const SymmetricMatrix<ElementType>& matrix, size_t column_index, size_t nonzero_index):
         matrix(matrix), column_index(column_index), nonzero_index(nonzero_index) { }

      [[nodiscard]] std::tuple<size_t, size_t, ElementType> operator*() const {
         return this->matrix.dereference_iterator(this->column_index, this->nonzero_index);
      }

      iterator& operator++() {
         this->matrix.increment_iterator(this->column_index, this->nonzero_index);
         return *this;
      }

      friend bool operator!=(const iterator& a, const iterator& b) {
         return &a.matrix != &b.matrix || a.column_index != b.column_index || a.nonzero_index != b.nonzero_index;
      }

   protected:
      const SymmetricMatrix<ElementType>& matrix;
      size_t column_index;
      size_t nonzero_index;
   };

   using value_type = ElementType;

   size_t dimension;
   size_t number_nonzeros{0};
   size_t capacity;

   SymmetricMatrix(size_t dimension, size_t original_capacity, bool use_regularization);
   virtual ~SymmetricMatrix() = default;

   virtual void reset();

   [[nodiscard]] size_t number_rows() const;
   [[nodiscard]] size_t number_columns() const;

   template <typename VectorType, typename ResultType>
   void product(const VectorType& vector, ResultType& result) const;

   [[nodiscard]] ElementType quadratic_product(const std::vector<ElementType>& x, const std::vector<ElementType>& y) const;

   // build the matrix incrementally
   virtual void insert(ElementType element, IndexType row_index, IndexType column_index) = 0;
   // this method will be used by the CSCSymmetricMatrix subclass
   virtual void finalize_column(IndexType column_index) = 0;
   [[nodiscard]] virtual ElementType smallest_diagonal_entry() const = 0;
   virtual void set_regularization(const std::function<ElementType(size_t /*index*/)>& regularization_function) = 0;

   [[nodiscard]] iterator begin() const { return iterator(*this, 0, 0); }
   [[nodiscard]] iterator end() const { return iterator(*this, this->dimension, this->number_nonzeros); }

   [[nodiscard]] const ElementType* data_raw_pointer() const;

   virtual void print(std::ostream& stream) const = 0;
   template <typename I, typename E>
   friend std::ostream& operator<<(std::ostream& stream, const SymmetricMatrix<I, E>& matrix);

protected:
   std::vector<ElementType> entries{};
   // regularization
   const bool use_regularization;

   // virtual iterator functions
   [[nodiscard]] virtual std::tuple<size_t, size_t, ElementType> dereference_iterator(size_t column_index, size_t nonzero_index) const = 0;
   virtual void increment_iterator(size_t& column_index, size_t& nonzero_index) const = 0;
};

// implementation

template <typename IndexType, typename ElementType>
SymmetricMatrix<IndexType, ElementType>::SymmetricMatrix(size_t dimension, size_t original_capacity, bool use_regularization) :
      dimension(dimension),
      // if regularization is used, allocate the necessary space
      capacity(original_capacity + (use_regularization ? dimension : 0)),
      use_regularization(use_regularization) {
   this->entries.reserve(this->capacity);
}

template <typename IndexType, typename ElementType>
void SymmetricMatrix<IndexType, ElementType>::reset() {
   this->number_nonzeros = 0;
   this->entries.clear();
}

template <typename IndexType, typename ElementType>
size_t SymmetricMatrix<IndexType, ElementType>::number_rows() const {
   return this->dimension;
}

template <typename IndexType, typename ElementType>
size_t SymmetricMatrix<IndexType, ElementType>::number_columns() const {
   return this->dimension;
}

template <typename IndexType, typename ElementType>
template <typename VectorType, typename ResultType>
void SymmetricMatrix<IndexType, ElementType>::product(const VectorType& vector, ResultType& result) const {
   this->for_each([&](IndexType row_index, IndexType column_index, double entry) {
      result[row_index] += entry*vector[column_index];
      // off-diagonal element on the other side of the diagonal
      if (row_index != column_index) {
         result[column_index] += entry*vector[row_index];
      }
   });
}

template <typename IndexType, typename ElementType>
ElementType SymmetricMatrix<IndexType, ElementType>::quadratic_product(const std::vector<ElementType>& x, const std::vector<ElementType>& y) const {
   assert(x.size() == y.size() && "SymmetricMatrix::quadratic_product: the two vectors x and y do not have the same size");

   ElementType result = ElementType(0);
   for (const auto [row_index, column_index, element]: *this) {
      result += (row_index == column_index ? ElementType(1) : ElementType(2)) * element * x[row_index] * y[column_index];
   }
   return result;
}

template <typename IndexType, typename ElementType>
const ElementType* SymmetricMatrix<IndexType, ElementType>::data_raw_pointer() const {
   return this->entries.data();
}

template <typename IndexType, typename ElementType>
std::ostream& operator<<(std::ostream& stream, const SymmetricMatrix<IndexType, ElementType>& matrix) {
   stream << "Dimension: " << matrix.dimension << ", number of nonzeros: " << matrix.number_nonzeros << '\n';
   matrix.print(stream);
   return stream;
}

#endif // UNO_SYMMETRICMATRIX_H
