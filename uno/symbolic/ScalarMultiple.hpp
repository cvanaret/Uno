// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SCALARMULTIPLE_H
#define UNO_SCALARMULTIPLE_H

namespace uno {
   // stores the expression (factor * expression) symbolically
   template <typename Expression>
   class ScalarMultiple {
   public:
      using value_type = typename std::remove_reference_t<Expression>::value_type;

      ScalarMultiple(value_type factor, Expression&& expression): factor(factor), expression(std::forward<Expression>(expression)) { }

      [[nodiscard]] constexpr size_t size() const { return this->expression.size(); }
      [[nodiscard]] value_type operator[](size_t index) const {
         return (this->factor == value_type(0)) ? value_type(0) : this->factor * this->expression[index];
      }

   protected:
      const value_type factor;
      Expression expression;
   };

   // free function
   template <typename Expression, typename ElementType = typename Expression::value_type,
         typename std::enable_if<std::is_arithmetic_v<ElementType>, int>::type = 0>
   inline ScalarMultiple<Expression> operator*(ElementType factor, Expression&& expression) {
      return ScalarMultiple<Expression>(factor, std::forward<Expression>(expression));
   }
} // namespace

#endif // UNO_SCALARMULTIPLE_H