// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ExactHessian.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "model/Model.hpp"
#include "tools/Logger.hpp"

namespace uno {
   size_t ExactHessian::compute_number_hessian_nonzeros(const Model& model) const {
      return model.number_hessian_nonzeros();
   }

   void ExactHessian::evaluate(const Model& model, const Vector<double>& primal_variables, double objective_multiplier,
         const Vector<double>& constraint_multipliers, SymmetricMatrix<size_t, double>& hessian) {
      DEBUG << "Evaluating model Hessian\n";
      model.evaluate_lagrangian_hessian(primal_variables, objective_multiplier, constraint_multipliers, hessian);
      this->evaluation_count++;
   }
} // namespace