// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INDICATOR_H
#define UNO_INDICATOR_H

// Indicator: componentwise evaluation of a predicate
template <typename E, typename Predicate>
class Indicator {
public:
   using value_type = typename std::remove_reference_t<E>::value_type;

   Indicator(E&& expression, double constant, Predicate predicate):
      expression(std::forward<E>(expression)), constant(constant), predicate(predicate) {
   }
   [[nodiscard]] size_t size() const { return this->expression.size(); }
   [[nodiscard]] typename Indicator::value_type operator[](size_t index) const {
      return this->predicate(this->expression[index], this->constant) ? 1. : 0.;
   }

protected:
   E expression;
   const double constant;
   Predicate* predicate; // function pointer
};

// free functions
inline bool less_than(double a, double b) { return a <= b; }
inline bool strictly_less_than(double a, double b) { return a < b; }
inline bool greater_than(double a, double b) { return a >= b; }
inline bool strictly_greater_than(double a, double b) { return a > b; }

template <typename E>
Indicator<E, decltype(less_than)> operator<=(E&& expression, double constant) {
   return {std::forward<E>(expression), constant, less_than};
}

template <typename E>
Indicator<E, decltype(strictly_less_than)> operator<(E&& expression, double constant) {
   return {std::forward<E>(expression), constant, strictly_less_than};
}

template <typename E>
Indicator<E, decltype(greater_than)> operator>=(E&& expression, double constant) {
   return {std::forward<E>(expression), constant, greater_than};
}

template <typename E>
Indicator<E, decltype(strictly_greater_than)> operator>(E&& expression, double constant) {
   return {std::forward<E>(expression), constant, strictly_greater_than};
}

#endif // UNO_INDICATOR_H