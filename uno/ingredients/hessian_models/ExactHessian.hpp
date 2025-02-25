// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_EXACTHESSIAN_H
#define UNO_EXACTHESSIAN_H

#include "HessianModel.hpp"

namespace uno {
   // exact Hessian
   class ExactHessian : public HessianModel {
   public:
      ExactHessian();

      void evaluate(const OptimizationProblem& problem, const Vector<double>& primal_variables, const Vector<double>& constraint_multipliers,
         SymmetricMatrix<size_t, double>& hessian) override;
   };
} // namespace

#endif // UNO_EXACTHESSIAN_H