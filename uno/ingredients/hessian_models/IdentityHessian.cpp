// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "IdentityHessian.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "model/Model.hpp"

namespace uno {
   void IdentityHessian::evaluate(const Model& model, const Vector<double>& /*primal_variables*/, double /*objective_multiplier*/,
         const Vector<double>& /*constraint_multipliers*/, SymmetricMatrix<size_t, double>& hessian) {
      hessian.reset();
      for (size_t variable_index: Range(model.number_variables)) {
         hessian.insert(1., variable_index, variable_index);
         hessian.finalize_column(variable_index);
      }
   }
} // namespace