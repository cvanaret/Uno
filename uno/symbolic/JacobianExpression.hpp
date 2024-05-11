// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_JACOBIANEXPRESSION_H
#define UNO_JACOBIANEXPRESSION_H

#include "linear_algebra/Matrix.hpp"
#include "optimization/Iterate.hpp"
#include "reformulation/OptimizationProblem.hpp"

class JacobianExpression {
public:
   using value_type = double;

   JacobianExpression(const OptimizationProblem& problem, Iterate& iterate): problem(problem), iterate(iterate) { }
   [[nodiscard]] size_t number_rows() const { return this->problem.number_constraints; }
   [[nodiscard]] size_t number_columns() const { return this->problem.number_variables; }
   void for_each(const std::function<void (size_t row_index, size_t column_index, double element)>& f) const;
   void evaluate(Matrix<double>& matrix);

protected:
   const OptimizationProblem& problem;
   Iterate& iterate;
};

inline void JacobianExpression::for_each(const std::function<void (size_t /*row_index*/, size_t /*column_index*/, double /*element*/)>& f) const {
   //RectangularMatrix<double> constraint_jacobian(this->number_rows(), this->number_columns());
   RectangularMatrix<double> constraint_jacobian(this->number_rows());
   for (auto& row: constraint_jacobian) {
      row.reserve(this->number_columns());
   }
   this->problem.evaluate_constraint_jacobian(this->iterate, constraint_jacobian);

   // go through the elements
   for (size_t constraint_index: Range(this->number_rows())) {
      constraint_jacobian[constraint_index].for_each([&](size_t variable_index, double derivative) {
         f(constraint_index, variable_index, derivative);
      });
   }
}

// free function
inline JacobianExpression jacobian(const OptimizationProblem& problem, Iterate& iterate) {
   return {problem, iterate};
}

inline void JacobianExpression::evaluate(Matrix<double>& /*matrix*/) {
   std::cout << "JacobianExpression::evaluate\n";
}

#endif // UNO_JACOBIANEXPRESSION_H