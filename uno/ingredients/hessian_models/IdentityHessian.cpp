// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "IdentityHessian.hpp"
#include "model/Model.hpp"
#include "tools/Logger.hpp"

namespace uno {
   bool IdentityHessian::has_hessian_operator(const Model& /*model*/) const {
      return true;
   }

   bool IdentityHessian::has_hessian_matrix(const Model& /*model*/) const {
      return true;
   }

   bool IdentityHessian::has_curvature(const Model& /*model*/) const {
      return true;
   }

   size_t IdentityHessian::number_nonzeros(const Model& model) const {
      return model.number_variables;
   }

   void IdentityHessian::compute_sparsity(const Model& model, int* row_indices, int* column_indices, int solver_indexing) const {
      // diagonal structure
      for (size_t variable_index: Range(model.number_variables)) {
         row_indices[variable_index] = static_cast<int>(variable_index) + solver_indexing;
         column_indices[variable_index] = static_cast<int>(variable_index) + solver_indexing;
      }
   }

   bool IdentityHessian::is_positive_definite() const {
      return true;
   }

   void IdentityHessian::initialize(const Model& /*model*/) {
   }

   void IdentityHessian::initialize_statistics(Statistics& /*statistics*/, const Options& /*options*/) const {
   }

   void IdentityHessian::notify_accepted_iterate(const Iterate& /*current_iterate*/, const Iterate& /*trial_iterate*/) {
      // do nothing
   }

   void IdentityHessian::evaluate_hessian(Statistics& /*statistics*/, const Model& model, const Vector<double>& /*primal_variables*/,
         double /*objective_multiplier*/, const Vector<double>& /*constraint_multipliers*/, double* hessian_values) {
      DEBUG << "Setting identity Hessian\n";
      for (size_t variable_index: Range(model.number_variables)) {
         hessian_values[variable_index] = 1.;
      }
   }

   void IdentityHessian::compute_hessian_vector_product(const Model& model, const double* /*x*/, const double* vector,
         double /*objective_multiplier*/, const Vector<double>& /*constraint_multipliers*/, double* result) {
      for (size_t variable_index: Range(model.number_variables)) {
         result[variable_index] = vector[variable_index];
      }
   }

   std::string IdentityHessian::get_name() const {
      return "identity";
   }
} // namespace