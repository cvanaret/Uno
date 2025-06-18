// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ZeroHessian.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "model/Model.hpp"

namespace uno {
   void ZeroHessian::initialize(const Model& /*model*/) {
   }

   size_t ZeroHessian::number_nonzeros(const Model& /*model*/) const {
      return 0;
   }

   bool ZeroHessian::is_positive_definite() const {
      return false;
   }

   void ZeroHessian::evaluate_hessian(Statistics& /*statistics*/, const Model& model, const Vector<double>& /*primal_variables*/,
         double /*objective_multiplier*/, const Vector<double>& /*constraint_multipliers*/, SymmetricMatrix<size_t, double>& hessian) {
      hessian.set_dimension(model.number_variables);
      hessian.reset();
      for (size_t variable_index: Range(model.number_variables)) {
         hessian.finalize_column(variable_index);
      }
   }

   void ZeroHessian::compute_hessian_vector_product(const Model& /*model*/, const Vector<double>& /*vector*/, double /*objective_multiplier*/,
         const Vector<double>& /*constraint_multipliers*/, Vector<double>& result) {
      result.fill(0.);
   }

   std::string ZeroHessian::get_name() const {
      return "zero";
   }
}