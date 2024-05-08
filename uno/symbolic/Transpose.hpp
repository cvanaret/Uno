// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_TRANSPOSE_H
#define UNO_TRANSPOSE_H

// Transpose: symbolic transpose of an expression
template <typename E>
class Transpose {
public:
   using value_type = typename std::remove_reference_t<E>::value_type;

   explicit Transpose(E&& expression): expression(std::forward<E>(expression)) {
   }
   [[nodiscard]] size_t size() const { return this->expression.size(); }
   [[nodiscard]] typename Transpose::value_type operator[](const std::pair<size_t, size_t>& indices) const {
      // switch the order of the indices
      return this->expression[{indices.second, indices.first}];
   }

protected:
   E expression;
};

template <typename E>
inline Transpose<E> transpose(E&& expression) {
   return {std::forward<E>(expression)};
}

#endif // UNO_TRANSPOSE_H