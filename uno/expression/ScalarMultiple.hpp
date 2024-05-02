// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SCALARMULTIPLE_H
#define UNO_SCALARMULTIPLE_H

// stores the expression (factor * expression) symbolically
template <typename ExpressionType>
class ScalarMultiple {
public:
   ScalarMultiple(double factor, ExpressionType&& expression): factor(factor), expression(std::forward<ExpressionType>(expression)) { }
   [[nodiscard]] size_t size() const { return this->expression.size(); }
   double operator[](size_t index) const { return (this->factor == 0.) ? 0. : this->factor * this->expression[index]; }

protected:
   const double factor;
   ExpressionType expression;
};

// free function
template <typename ExpressionType>
inline ScalarMultiple<ExpressionType> operator*(double factor, ExpressionType&& expression) {
   return ScalarMultiple<ExpressionType>(factor, std::forward<ExpressionType>(expression));
}

#endif // UNO_SCALARMULTIPLE_H