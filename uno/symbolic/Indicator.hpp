// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INDICATOR_H
#define UNO_INDICATOR_H

// Indicator: componentwise evaluation of a predicate
template <typename Expression, typename Predicate>
class Indicator {
public:
   using value_type = typename std::remove_reference_t<Expression>::value_type;

   Indicator(Expression&& expression, value_type constant, Predicate predicate):
      expression(std::forward<Expression>(expression)), constant(constant), predicate(predicate) {
   }
   [[nodiscard]] size_t size() const { return this->expression.size(); }
   [[nodiscard]] value_type operator[](size_t index) const {
      return this->predicate(this->expression[index], this->constant) ? value_type(1) : value_type(0);
   }

protected:
   Expression expression;
   const value_type constant;
   Predicate* predicate; // function pointer
};

// free functions
template <typename ElementType>
inline bool less_than(ElementType a, ElementType b) { return a <= b; }
template <typename ElementType>
inline bool strictly_less_than(ElementType a, ElementType b) { return a < b; }
template <typename ElementType>
inline bool greater_than(ElementType a, ElementType b) { return a >= b; }
template <typename ElementType>
inline bool strictly_greater_than(ElementType a, ElementType b) { return a > b; }

template <typename Expression, typename ElementType = typename Expression::value_type>
Indicator<Expression, decltype(less_than<ElementType>)> operator<=(Expression&& expression, ElementType constant) {
   return {std::forward<Expression>(expression), constant, less_than<ElementType>};
}

template <typename Expression, typename ElementType = typename Expression::value_type>
Indicator<Expression, decltype(strictly_less_than<ElementType>)> operator<(Expression&& expression, ElementType constant) {
   return {std::forward<Expression>(expression), constant, strictly_less_than<ElementType>};
}

template <typename Expression, typename ElementType = typename Expression::value_type>
Indicator<Expression, decltype(greater_than<ElementType>)> operator>=(Expression&& expression, ElementType constant) {
   return {std::forward<Expression>(expression), constant, greater_than<ElementType>};
}

template <typename Expression, typename ElementType = typename Expression::value_type>
Indicator<Expression, decltype(strictly_greater_than<ElementType>)> operator>(Expression&& expression, ElementType constant) {
   return {std::forward<Expression>(expression), constant, strictly_greater_than<ElementType>};
}

#endif // UNO_INDICATOR_H