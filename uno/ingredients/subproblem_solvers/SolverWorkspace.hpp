// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_EVALUATIONSPACE_H
#define UNO_EVALUATIONSPACE_H

#include "../../optimization/Evaluations.hpp"

namespace uno {
   // forward declaration
   template <typename ElementType>
   class DirectSymmetricIndefiniteLinearSolver;
   class Model;
   class OptimizationProblem;
   class Statistics;
   class Subproblem;
   template <typename ElementType>
   class Vector;
   class WarmstartInformation;

   class SolverWorkspace {
   public:
      SolverWorkspace() = default;
      virtual ~SolverWorkspace() = default;

      [[nodiscard]] virtual double compute_hessian_quadratic_product(const Subproblem& subproblem, const Vector<double>& vector) const = 0;
   };
} // namespace

#endif // UNO_EVALUATIONSPACE_H