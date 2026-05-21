// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_UNARYNEGATION_H
#define UNO_UNARYNEGATION_H

#include <cstddef>
#include "symbolic_traits.hpp"

namespace uno {
   // stores the expression -expression symbolically
   // limited to types that possess value_type
   // https://stackoverflow.com/questions/11055923/stdenable-if-parameter-vs-template-parameter
   template <typename Expression>
   class UnaryNegation {
   public:
      using value_type = typename std::remove_reference_t<Expression>::value_type;

      explicit UnaryNegation(Expression&& expression): expression(std::forward<Expression>(expression)) { }

      [[nodiscard]] constexpr size_t size() const noexcept {
         return this->expression.size();
      }

      [[nodiscard]] constexpr value_type operator[](size_t index) const noexcept {
         return -this->expression[index];
      }

   protected:
      storage_t<Expression> expression;
   };

   // free function
   template <typename Expression>
   UnaryNegation<Expression> operator-(Expression&& expression) {
      return UnaryNegation<Expression>(std::forward<Expression>(expression));
   }
} // namespace

#endif // UNO_UNARYNEGATION_H
