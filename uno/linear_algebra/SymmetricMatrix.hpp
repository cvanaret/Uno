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
   virtual void set_regularization(const std::function<ElementType(IndexType /*index*/)>& regularization_function) = 0;

   [[nodiscard]] const ElementType* data_raw_pointer() const;

   virtual void print(std::ostream& stream) const = 0;
   template <typename I, typename E>
   friend std::ostream& operator<<(std::ostream& stream, const SymmetricMatrix<I, E>& matrix);

protected:
   std::vector<ElementType> entries{};
   // regularization
   const bool use_regularization;
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
   this->for_each([&](IndexType row_index, IndexType column_index, ElementType entry) {
      result += (row_index == column_index ? ElementType(1) : ElementType(2)) * entry * x[row_index] * y[column_index];
   });
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
