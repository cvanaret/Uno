// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ZeroHessian.hpp"
#include "symbolic/Range.hpp"

namespace uno {
   ZeroHessian::ZeroHessian(size_t number_variables): HessianModel("zero"), number_variables(number_variables) {
   }

   bool ZeroHessian::has_hessian_operator() const {
      return true;
   }

   bool ZeroHessian::has_hessian_matrix() const {
      return true;
   }

   bool ZeroHessian::has_curvature() const {
      return false;
   }

   size_t ZeroHessian::number_nonzeros() const {
      return 0;
   }

   void ZeroHessian::compute_sparsity(uno_int* /*row_indices*/, uno_int* /*column_indices*/, uno_int /*solver_indexing*/) const {
      // empty structure
   }

   bool ZeroHessian::is_positive_definite() const {
      return false;
   }

   void ZeroHessian::evaluate_hessian(Statistics& /*statistics*/, const Vector<double>& /*primal_variables*/,
         double /*objective_multiplier*/, const Vector<double>& /*constraint_multipliers*/, double* /*hessian_values*/) {
      // do nothing
   }

   void ZeroHessian::compute_hessian_vector_product(const double* /*x*/, const double* /*vector*/,
         double /*objective_multiplier*/, const Vector<double>& /*constraint_multipliers*/, double* result) {
      for (size_t variable_index: Range(this->number_variables)) {
         result[variable_index] = 0.;
      }
   }
}