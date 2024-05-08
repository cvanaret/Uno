// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_TRANSPOSE_H
#define UNO_TRANSPOSE_H

#include "MatrixExpression.hpp"

// Transpose: symbolic transpose of an expression
template <typename E>
class Transpose {
public:
   using value_type = typename std::remove_reference_t<E>::value_type;

   explicit Transpose(E&& expression): expression(std::forward<E>(expression)) { }

   [[nodiscard]] size_t size() const { return this->expression.size(); }

   [[nodiscard]] typename Transpose::value_type operator[](const std::pair<size_t, size_t>& indices) const {
      // switch the order of the indices
      return this->expression[{indices.second, indices.first}];
   }

   template <typename VectorType, typename ResultType>
   void product(const VectorType& /*vector*/, ResultType& /*result*/) const {
      throw std::runtime_error("Transpose::product is not implemented yet");
   }

protected:
   E expression;
};

// free function
template <typename E>
inline auto transpose(E&& expression) {
   // any diagonal matrix is its own transpose
   if constexpr (std::is_base_of_v<DiagonalMatrix, std::remove_reference_t<E>>) {
      return std::forward<E>(expression);
   }
   else {
      return Transpose<E>(std::forward<E>(expression));
   }
}

#endif // UNO_TRANSPOSE_H