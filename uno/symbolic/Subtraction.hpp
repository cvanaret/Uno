// Copyright (c) 2025-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SUBTRACTION_H
#define UNO_SUBTRACTION_H

#include <cstddef>
#include "symbolic_traits.hpp"

namespace uno {
   // stores the expression (left - right) symbolically
   // limited to types that possess value_type
   // https://stackoverflow.com/questions/11055923/stdenable-if-parameter-vs-template-parameter
   template <typename L, typename R,
      std::enable_if_t<std::is_same_v<typename std::remove_reference_t<L>::value_type,
                                      typename std::remove_reference_t<R>::value_type>, int> = 0>
   class Subtraction {
   public:
      using value_type = typename std::remove_reference_t<L>::value_type;

      Subtraction(L&& left, R&& right): left(std::forward<L>(left)), right(std::forward<R>(right)) { }

      [[nodiscard]] constexpr size_t size() const noexcept {
         return this->left.size();
      }

      [[nodiscard]] constexpr value_type operator[](size_t index) const noexcept {
         return this->left[index] - this->right[index];
      }

      UNO_FORWARD_ACCESSOR(get_left, this->left)

      UNO_FORWARD_ACCESSOR(get_right, this->right)

   protected:
      storage_t<L> left;
      storage_t<R> right;
   };

   // free function
   template <typename E1, typename E2>
   Subtraction<E1, E2> operator-(E1&& left, E2&& right) {
      return {std::forward<E1>(left), std::forward<E2>(right)};
   }
} // namespace

#endif // UNO_SUBTRACTION_H