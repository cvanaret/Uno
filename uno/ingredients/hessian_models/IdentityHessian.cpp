// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "IdentityHessian.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "options/Options.hpp"

namespace uno {
   void IdentityHessian::evaluate(const OptimizationProblem& problem, const Vector<double>& /*primal_variables*/,
         const Vector<double>& /*constraint_multipliers*/, SymmetricMatrix<size_t, double>& hessian) {
      hessian.set_dimension(problem.number_variables);
      for (size_t variable_index: Range(problem.number_variables)) {
         hessian.insert(1., variable_index, variable_index);
      }
   }
} // namespace