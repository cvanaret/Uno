// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "IdentityHessian.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "model/Model.hpp"
#include "tools/Logger.hpp"

namespace uno {
   void IdentityHessian::initialize(const Model& /*model*/) {
      // do nothing
   }

   size_t IdentityHessian::number_nonzeros(const Model& model) const {
      return model.number_variables;
   }

   void IdentityHessian::evaluate_hessian(Statistics& /*statistics*/, const Model& model, const Vector<double>& /*primal_variables*/,
         double /*objective_multiplier*/, const Vector<double>& /*constraint_multipliers*/, SymmetricMatrix<size_t, double>& hessian) {
      DEBUG << "Setting identity Hessian\n";
      hessian.reset();
      for (size_t variable_index: Range(model.number_variables)) {
         hessian.insert(1., variable_index, variable_index);
         hessian.finalize_column(variable_index);
      }
   }

   void IdentityHessian::compute_hessian_vector_product(const Model& model, const Vector<double>& vector, double /*objective_multiplier*/,
         const Vector<double>& /*constraint_multipliers*/, Vector<double>& result) {
      result.fill(0.);
      for (size_t variable_index: Range(model.number_variables)) {
         result[variable_index] = vector[variable_index];
      }
   }
} // namespace