// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MATRIXEXPRESSION_H
#define UNO_MATRIXEXPRESSION_H

#include "Expression.hpp"
#include "tools/Range.hpp"

class IdentityMatrix: public Expression {
public:
   using value_type = double;

   explicit IdentityMatrix(size_t dimension): Expression(), dimension(dimension) {
   }
   [[nodiscard]] size_t number_rows() const { return this->dimension; }
   [[nodiscard]] size_t number_columns() const { return this->dimension; }
   void for_each(const std::function<void (size_t row_index, size_t column_index, double element)>& f) const override {
      for (size_t index: Range(this->dimension)) {
         f(index, index, 1.);
      }
   }
   template <typename VectorType, typename ResultType>
   void product(const VectorType& vector, ResultType& result) const {
      for (size_t index: Range(result.size())) {
         result[index] += vector[index];
      }
   }

protected:
   const size_t dimension;
};

// free function
inline IdentityMatrix identity(size_t dimension) {
   return IdentityMatrix(dimension);
}

class DiagonalMatrix: public Expression {
public:
   using value_type = double;

   DiagonalMatrix(size_t dimension, double factor): Expression(), dimension(dimension), factor(factor) {
   }
   [[nodiscard]] size_t number_rows() const { return this->dimension; }
   [[nodiscard]] size_t number_columns() const { return this->dimension; }
   void for_each(const std::function<void (size_t row_index, size_t column_index, double element)>& f) const override {
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
   const double factor;
};

// free functions
inline DiagonalMatrix zeros(size_t dimension) {
   return DiagonalMatrix(dimension, 0.);
}

inline DiagonalMatrix diagonal(size_t dimension, double factor) {
   return DiagonalMatrix(dimension, factor);
}

#endif // UNO_MATRIXEXPRESSION_H