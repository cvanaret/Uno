// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MATRIXEXPRESSION_H
#define UNO_MATRIXEXPRESSION_H

#include "Range.hpp"

// square diagonal matrix
template <typename ElementType>
class DiagonalMatrix {
public:
   using value_type = ElementType;

   DiagonalMatrix(size_t dimension, ElementType factor): dimension(dimension), factor(factor) { }
   [[nodiscard]] size_t number_rows() const { return this->dimension; }
   [[nodiscard]] size_t number_columns() const { return this->dimension; }
   void for_each(const std::function<void (size_t row_index, size_t column_index, ElementType element)>& f) const {
      for (size_t index: Range(this->dimension)) {
         f(index, index, this->factor);
      }
   }
   template <typename VectorType, typename ResultType>
   void product(const VectorType& vector, ResultType& result) const {
      for (size_t index: Range(result.size())) {
         result[index] += this->factor * vector[index];
      }
   }

protected:
   const size_t dimension;
   const ElementType factor;
};

// free function
template <typename ElementType>
inline DiagonalMatrix<ElementType> diagonal(size_t dimension, ElementType factor) {
   return {dimension, factor};
}

// identity matrix
template <typename ElementType>
class IdentityMatrix: public DiagonalMatrix<ElementType> {
public:
   explicit IdentityMatrix(size_t dimension): DiagonalMatrix<ElementType>(dimension, ElementType(1)) { }
};

// free function
template <typename ElementType>
inline IdentityMatrix<ElementType> identity(size_t dimension) {
   return IdentityMatrix<ElementType>(dimension);
}

// zero matrix: maintains sparsity
template <typename ElementType>
class ZeroMatrix {
public:
   using value_type = ElementType;

   ZeroMatrix(size_t number_rows, size_t number_columns): number_of_rows(number_rows), number_of_columns(number_columns) { }
   [[nodiscard]] size_t number_rows() const { return this->number_of_rows; }
   [[nodiscard]] size_t number_columns() const { return this->number_of_columns; }
   void for_each(const std::function<void (size_t, size_t, ElementType)>& /*f*/) const { } // do nothing
   template <typename VectorType, typename ResultType>
   void product(const VectorType& /*vector*/, ResultType& /*result*/) const { } // do nothing

protected:
   const size_t number_of_rows, number_of_columns;
};

// free functions
template <typename ElementType>
inline ZeroMatrix<ElementType> zeros(size_t dimension) {
   return {dimension, dimension};
}

template <typename ElementType>
inline ZeroMatrix<ElementType> zeros(size_t number_rows, size_t number_columns) {
   return {number_rows, number_columns};
}

#endif // UNO_MATRIXEXPRESSION_H
