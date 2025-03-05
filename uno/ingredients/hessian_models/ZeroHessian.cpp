// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ZeroHessian.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "model/Model.hpp"

namespace uno {
   void ZeroHessian::evaluate(const Model& model, const Vector<double>& /*primal_variables*/, double /*objective_multiplier*/,
         const Vector<double>& /*constraint_multipliers*/, SymmetricMatrix<size_t, double>& hessian) {
      hessian.set_dimension(model.number_variables);
      hessian.reset();
   }
}