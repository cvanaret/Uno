// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MULTIPLICATION_H
#define UNO_MULTIPLICATION_H

#include "symbolic_traits.hpp"

namespace uno {
   // stores the expression (left * right) symbolically
   template <typename L, typename R>
   class Multiplication {
   public:
      using value_type = typename std::remove_reference_t<L>::value_type;

      Multiplication(L&& left, R&& right): left(std::forward<L>(left)), right(std::forward<R>(right)) { }

      [[nodiscard]] constexpr decltype(auto) get_left() const noexcept {
         return this->left;
      }

      [[nodiscard]] constexpr decltype(auto) get_right() const noexcept {
         return this->right;
      }

   protected:
      storage_t<L> left;
      storage_t<R> right;
   };

   // free function
   template <typename Matrix1, typename Matrix2,
         std::enable_if_t<std::is_same_v<std::remove_const_t<typename std::remove_reference_t<Matrix1>::value_type>,
                                         std::remove_const_t<typename std::remove_reference_t<Matrix2>::value_type>>, int> = 0,
      // Matrix1 and Matrix2 are both not arithmetic types
      std::enable_if_t<!std::is_arithmetic_v<std::remove_reference_t<Matrix1>>, int> = 0,
      std::enable_if_t<!std::is_arithmetic_v<std::remove_reference_t<Matrix2>>, int> = 0>
   Multiplication<Matrix1, Matrix2> operator*(Matrix1&& left, Matrix2&& right) {
      return {std::forward<Matrix1>(left), std::forward<Matrix2>(right)};
   }
} // namespace

#endif // UNO_MULTIPLICATION_H