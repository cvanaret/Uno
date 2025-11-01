// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "IdentityHessian.hpp"
#include "symbolic/Range.hpp"
#include "tools/Logger.hpp"

namespace uno {
   IdentityHessian::IdentityHessian(size_t number_variables): HessianModel("identity"), number_variables(number_variables) {
   }

   bool IdentityHessian::has_hessian_operator() const {
      return true;
   }

   bool IdentityHessian::has_hessian_matrix() const {
      return true;
   }

   bool IdentityHessian::has_curvature() const {
      return true;
   }

   size_t IdentityHessian::number_nonzeros() const {
      return this->number_variables;
   }

   void IdentityHessian::compute_sparsity(int* row_indices, int* column_indices, int solver_indexing) const {
      // diagonal structure
      for (size_t variable_index: Range(this->number_variables)) {
         row_indices[variable_index] = static_cast<int>(variable_index) + solver_indexing;
         column_indices[variable_index] = static_cast<int>(variable_index) + solver_indexing;
      }
   }

   bool IdentityHessian::is_positive_definite() const {
      return true;
   }

   void IdentityHessian::initialize_statistics(Statistics& /*statistics*/, const Options& /*options*/) const {
   }

   void IdentityHessian::notify_accepted_iterate(Iterate& /*current_iterate*/, Iterate& /*trial_iterate*/) {
   }

   void IdentityHessian::evaluate_hessian(Statistics& /*statistics*/, const Vector<double>& /*primal_variables*/,
         double /*objective_multiplier*/, const Vector<double>& /*constraint_multipliers*/, double* hessian_values) {
      DEBUG << "Setting identity Hessian\n";
      for (size_t variable_index: Range(this->number_variables)) {
         hessian_values[variable_index] = 1.;
      }
   }

   void IdentityHessian::compute_hessian_vector_product(const double* /*x*/, const double* vector,
         double /*objective_multiplier*/, const Vector<double>& /*constraint_multipliers*/, double* result) {
      for (size_t variable_index: Range(this->number_variables)) {
         result[variable_index] = vector[variable_index];
      }
   }
} // namespace