// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SUBTRACTION_H
#define UNO_SUBTRACTION_H

namespace uno {
   // stores the expression (expression1 - expression2) symbolically
   // limited to types that possess value_type
   // https://stackoverflow.com/questions/11055923/stdenable-if-parameter-vs-template-parameter
   template <typename E1, typename E2,
      typename std::enable_if_t<std::is_same_v<typename std::remove_reference_t<E1>::value_type, typename std::remove_reference_t<E2>::value_type>, int> = 0>
   class Subtraction {
   public:
      using value_type = typename std::remove_reference_t<E1>::value_type;

      Subtraction(E1&& expression1, E2&& expression2): expression1(std::forward<E1>(expression1)), expression2(std::forward<E2>(expression2)) { }

      [[nodiscard]] constexpr size_t size() const {
         return this->expression1.size();
      }

      [[nodiscard]] typename Subtraction::value_type operator[](size_t index) const {
         return this->expression1[index] - this->expression2[index];
      }

   protected:
      const E1 expression1;
      const E2 expression2;
   };

   // free function
   template <typename E1, typename E2>
   inline Subtraction<E1, E2> operator-(E1&& expression1, E2&& expression2) {
      return {std::forward<E1>(expression1), std::forward<E2>(expression2)};
   }
} // namespace

#endif // UNO_SUBTRACTION_H
