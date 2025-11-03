// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ExactHessian.hpp"
#include "model/Model.hpp"

namespace uno {
   ExactHessian::ExactHessian(const Model& model): HessianModel("exact"), model(model) {
   }

   bool ExactHessian::has_hessian_operator() const {
      return this->model.has_hessian_operator();
   }

   bool ExactHessian::has_hessian_matrix() const {
      return this->model.has_hessian_matrix();
   }

   bool ExactHessian::has_curvature() const {
      return (this->model.get_problem_type() != ProblemType::LINEAR);
   }

   size_t ExactHessian::number_nonzeros() const {
      return this->model.number_hessian_nonzeros();
   }

   void ExactHessian::compute_sparsity(uno_int* row_indices, uno_int* column_indices, uno_int solver_indexing) const {
      // Hessian sparsity of the model
      this->model.compute_hessian_sparsity(row_indices, column_indices, solver_indexing);
   }

   bool ExactHessian::is_positive_definite() const {
      return false;
   }

   void ExactHessian::evaluate_hessian(Statistics& /*statistics*/, const Vector<double>& primal_variables,
         double objective_multiplier, const Vector<double>& constraint_multipliers, double* hessian_values) {
      this->model.evaluate_lagrangian_hessian(primal_variables, objective_multiplier, constraint_multipliers, hessian_values);
      ++this->evaluation_count;
   }

   void ExactHessian::compute_hessian_vector_product(const double* x, const double* vector,
         double objective_multiplier, const Vector<double>& constraint_multipliers, double* result) {
      this->model.compute_hessian_vector_product(x, vector, objective_multiplier, constraint_multipliers, result);
      ++this->evaluation_count;
   }
} // namespace