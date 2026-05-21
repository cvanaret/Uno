// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SUM_H
#define UNO_SUM_H

#include "symbolic_traits.hpp"

namespace uno {
   // stores the expression (left + right) symbolically
   // limited to types that possess value_type
   // https://stackoverflow.com/questions/11055923/stdenable-if-parameter-vs-template-parameter
   template <typename E1, typename E2,
      std::enable_if_t<std::is_same_v<typename std::remove_reference_t<E1>::value_type,
                                      typename std::remove_reference_t<E2>::value_type>, int> = 0>
   class Sum {
   public:
      using value_type = typename std::remove_reference_t<E1>::value_type;

      Sum(E1&& left, E2&& right): left(std::forward<E1>(left)), right(std::forward<E2>(right)) { }

      [[nodiscard]] constexpr size_t size() const noexcept {
         return this->left.size();
      }

      [[nodiscard]] constexpr value_type operator[](size_t index) const noexcept {
         return this->left[index] + this->right[index];
      }

      [[nodiscard]] constexpr decltype(auto) get_left() const noexcept {
         return this->left;
      }

      [[nodiscard]] constexpr decltype(auto) get_right() const noexcept {
         return this->right;
      }

   protected:
      storage_t<E1> left;
      storage_t<E2> right;
   };

   // free function
   template <typename L, typename R>
   Sum<L, R> operator+(L&& left, R&& right) {
      return {std::forward<L>(left), std::forward<R>(right)};
   }
} // namespace

#endif // UNO_SUM_H