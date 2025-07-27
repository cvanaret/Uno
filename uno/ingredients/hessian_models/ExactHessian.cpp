// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ExactHessian.hpp"
#include "model/Model.hpp"

namespace uno {
   bool ExactHessian::has_implicit_representation() const {
      // As long as we use the ASL library ("solvers"), we need to form the explicit Hessian
      // The reason is that the ASL Hessian representation changes as soon as trial
      // iterates are evaluated. The variant "solvers2" should address the issue.
      return false;
   }

   bool ExactHessian::has_explicit_representation() const {
      return true;
   }

   bool ExactHessian::has_curvature(const Model& model) const {
      return (0 < model.number_hessian_nonzeros());
   }

   size_t ExactHessian::number_nonzeros(const Model& model) const {
      return model.number_hessian_nonzeros();
   }

   void ExactHessian::compute_structure(const Model& model, size_t* row_indices, size_t* column_indices, size_t solver_indexing) const {
      // Hessian structure of the model
      model.compute_hessian_structure(row_indices, column_indices, solver_indexing);
   }

   bool ExactHessian::is_positive_definite() const {
      return false;
   }

   void ExactHessian::initialize(const Model& /*model*/) {
   }

   void ExactHessian::evaluate_hessian(Statistics& /*statistics*/, const Model& model, const Vector<double>& primal_variables,
         double objective_multiplier, const Vector<double>& constraint_multipliers, Vector<double>& hessian_values) {
      model.evaluate_lagrangian_hessian(primal_variables, objective_multiplier, constraint_multipliers, hessian_values);
      this->evaluation_count++;
   }

   void ExactHessian::compute_hessian_vector_product(const Model& model, const double* vector, double objective_multiplier,
         const Vector<double>& constraint_multipliers, double* result) {
      model.compute_hessian_vector_product(vector, objective_multiplier, constraint_multipliers, result);
      this->evaluation_count++;
   }

   std::string ExactHessian::get_name() const {
      return "exact";
   }
} // namespace