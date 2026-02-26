// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MULTIPLICATION_H
#define UNO_MULTIPLICATION_H

#include <utility>

namespace uno {
   // stores the expression (expression1 * expression2) symbolically
   // limited to types that possess value_type
   // https://stackoverflow.com/questions/11055923/stdenable-if-parameter-vs-template-parameter
   template <typename E1, typename E2,
      typename std::enable_if_t<std::is_same_v<typename std::remove_reference_t<E1>::value_type, typename std::remove_reference_t<E2>::value_type>, int> = 0>
   class Multiplication {
   public:
      using value_type = typename std::remove_reference_t<E1>::value_type;

      Multiplication(E1&& expression1, E2&& expression2): expression1(std::forward<E1>(expression1)), expression2(std::forward<E2>(expression2)) { }

      [[nodiscard]] constexpr size_t size() const {
         return this->expression1.size();
      }

      [[nodiscard]] const E1& get_left() const {
         return this->expression1;
      }

      [[nodiscard]] const E2& get_right() const {
         return this->expression2;
      }

   protected:
      const E1 expression1;
      const E2 expression2;
   };

   // free function
   template <typename E1, typename E2>
   inline Multiplication<E1, E2> operator*(E1&& expression1, E2&& expression2) {
      return {std::forward<E1>(expression1), std::forward<E2>(expression2)};
   }
} // namespace

#endif // UNO_MULTIPLICATION_H