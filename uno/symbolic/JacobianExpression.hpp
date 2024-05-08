// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_JACOBIANEXPRESSION_H
#define UNO_JACOBIANEXPRESSION_H

#include "optimization/Iterate.hpp"
#include "reformulation/OptimizationProblem.hpp"

class JacobianExpression {
public:
   using value_type = double;

   JacobianExpression(const OptimizationProblem& problem, const Iterate& iterate): problem(problem), iterate(iterate) { }
   [[nodiscard]] size_t number_rows() const { return this->problem.number_variables; }
   [[nodiscard]] size_t number_columns() const { return this->problem.number_constraints; }
   void for_each(const std::function<void (size_t row_index, size_t column_index, double element)>& f) const;

protected:
   const OptimizationProblem& problem;
   const Iterate& iterate;
};

inline void JacobianExpression::for_each(const std::function<void (size_t /*row_index*/, size_t /*column_index*/, double /*element*/)>& /*f*/) const {
   /*this->problem.evaluate_constraint_jacobian(this->iterate, RectangularMatrix<double>);

   // go through the elements
   this->hessian_problem.hessian->for_each(f);*/
}

// free function
inline JacobianExpression jac(const OptimizationProblem& problem, const Iterate& iterate) {
   return {problem, iterate};
}

#endif // UNO_JACOBIANEXPRESSION_H