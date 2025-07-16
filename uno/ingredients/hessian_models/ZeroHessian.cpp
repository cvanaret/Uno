// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ZeroHessian.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "linear_algebra/Vector.hpp"
#include "model/Model.hpp"

namespace uno {
   bool ZeroHessian::has_implicit_representation() const {
      return true;
   }

   bool ZeroHessian::has_explicit_representation() const {
      return true;
   }

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
      hessian.reset();
      for (size_t variable_index: Range(model.number_variables)) {
         hessian.finalize_column(variable_index);
      }
   }

   void ZeroHessian::compute_hessian_vector_product(const Model& model, const double* /*vector*/, double /*objective_multiplier*/,
         const Vector<double>& /*constraint_multipliers*/, double* result) {
      for (size_t variable_index: Range(model.number_variables)) {
         result[variable_index] = 0.;
      }
   }

   std::string ZeroHessian::get_name() const {
      return "zero";
   }
}