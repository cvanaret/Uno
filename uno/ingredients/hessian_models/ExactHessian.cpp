// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ExactHessian.hpp"
#include "model/Model.hpp"

namespace uno {
   bool ExactHessian::has_implicit_representation(const Model& model) const {
      return model.has_implicit_hessian_representation();
   }

   bool ExactHessian::has_explicit_representation(const Model& model) const {
      return model.has_explicit_hessian_representation();
   }

   bool ExactHessian::has_curvature(const Model& model) const {
      return (0 < model.number_hessian_nonzeros());
   }

   size_t ExactHessian::number_nonzeros(const Model& model) const {
      return model.number_hessian_nonzeros();
   }

   void ExactHessian::compute_sparsity(const Model& model, int* row_indices, int* column_indices, int solver_indexing) const {
      // Hessian sparsity of the model
      model.compute_hessian_sparsity(row_indices, column_indices, solver_indexing);
   }

   bool ExactHessian::is_positive_definite() const {
      return false;
   }

   void ExactHessian::initialize(const Model& /*model*/) {
   }

   void ExactHessian::evaluate_hessian(Statistics& /*statistics*/, const Model& model, const Vector<double>& primal_variables,
         double objective_multiplier, const Vector<double>& constraint_multipliers, double* hessian_values) {
      model.evaluate_lagrangian_hessian(primal_variables, objective_multiplier, constraint_multipliers, hessian_values);
      ++this->evaluation_count;
   }

   void ExactHessian::compute_hessian_vector_product(const Model& model, const double* x, const double* vector,
         double objective_multiplier, const Vector<double>& constraint_multipliers, double* result) {
      model.compute_hessian_vector_product(x, vector, objective_multiplier, constraint_multipliers, result);
      ++this->evaluation_count;
   }

   std::string ExactHessian::get_name() const {
      return "exact";
   }
} // namespace