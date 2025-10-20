// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_EVALUATIONSPACE_H
#define UNO_EVALUATIONSPACE_H

namespace uno {
   // forward declaration
   template <typename ElementType>
   class DirectSymmetricIndefiniteLinearSolver;
   class Iterate;
   class OptimizationProblem;
   class Statistics;
   class Subproblem;
   template <typename ElementType>
   class Vector;
   class WarmstartInformation;

   class EvaluationSpace {
   public:
      EvaluationSpace() = default;
      virtual ~EvaluationSpace() = default;

      virtual void evaluate_constraint_jacobian(const OptimizationProblem& problem, Iterate& iterate) = 0;
      virtual void compute_constraint_jacobian_vector_product(const Vector<double>& vector, Vector<double>& result) const = 0;
      virtual void compute_constraint_jacobian_transposed_vector_product(const Vector<double>& vector,
         Vector<double>& result) const = 0;
      [[nodiscard]] virtual double compute_hessian_quadratic_product(const Subproblem& subproblem, const Vector<double>& vector) const = 0;
   };
} // namespace

#endif // UNO_EVALUATIONSPACE_H