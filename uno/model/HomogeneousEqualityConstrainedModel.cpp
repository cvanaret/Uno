// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "HomogeneousEqualityConstrainedModel.hpp"
#include "linear_algebra/VectorView.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/Range.hpp"

namespace uno {
   // Transform the problem into an equality-constrained problem with constraints c(x) = 0. This implies:
   // - inequality constraints get a slack
   // - equality constraints are shifted by their RHS
   HomogeneousEqualityConstrainedModel::HomogeneousEqualityConstrainedModel(const Model& original_model):
         Model(original_model.name + " -> equality constrained", original_model.number_variables +
            original_model.get_inequality_constraints().size(), original_model.number_constraints,
            original_model.optimization_sense, original_model.lagrangian_sign_convention, original_model.base_indexing),
         model(original_model),
         // all constraints are equality constraints
         equality_constraints(Range(this->number_constraints)),
         inequality_constraints(Range(0)),
         slacks(this->model.get_inequality_constraints().size()),
         variables_lower_bounds(this->number_variables),
         variables_upper_bounds(this->number_variables),
         constraints_lower_bounds(this->number_constraints, 0.),
         constraints_upper_bounds(this->number_constraints, 0.) {
      // copy the original variables bounds
      view(this->variables_lower_bounds, 0, original_model.number_variables) = original_model.get_variables_lower_bounds();
      view(this->variables_upper_bounds, 0, original_model.number_variables) = original_model.get_variables_upper_bounds();
      // register the inequality constraint of each slack
      size_t inequality_index = 0;
      for (const size_t constraint_index: this->model.get_inequality_constraints()) {
         const size_t slack_variable_index = this->model.number_variables + inequality_index;
         this->slacks.insert(constraint_index, slack_variable_index);
         // bounds of the slack
         this->variables_lower_bounds[slack_variable_index] = this->model.get_constraints_lower_bounds()[constraint_index];
         this->variables_upper_bounds[slack_variable_index] = this->model.get_constraints_upper_bounds()[constraint_index];
         ++inequality_index;
      }
   }

   ProblemType HomogeneousEqualityConstrainedModel::get_problem_type() const {
      return this->model.get_problem_type();
   }

   bool HomogeneousEqualityConstrainedModel::has_jacobian_operator() const {
      return this->model.has_jacobian_operator();
   }

   bool HomogeneousEqualityConstrainedModel::has_jacobian_transposed_operator() const {
      return this->model.has_jacobian_transposed_operator();
   }

   bool HomogeneousEqualityConstrainedModel::has_hessian_operator() const {
      return this->model.has_hessian_operator();
   }
   
   bool HomogeneousEqualityConstrainedModel::has_hessian_matrix() const {
      return this->model.has_hessian_matrix();
   }

   double HomogeneousEqualityConstrainedModel::evaluate_objective(const Vector<double>& x) const {
      return this->model.evaluate_objective(x);
   }

   void HomogeneousEqualityConstrainedModel::evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const {
      this->model.evaluate_constraints(x, constraints);
      // inequality constraints: add the slacks
      for (const auto [constraint_index, slack_index]: this->get_slacks()) {
         constraints[constraint_index] -= x[slack_index];
      }

      // equality constraints: make sure they are homogeneous (c(x) = 0)
      for (const size_t constraint_index: this->model.get_equality_constraints()) {
         const double fixed_bound = this->model.get_constraints_lower_bounds()[constraint_index];
         constraints[constraint_index] -= fixed_bound;
      }
   }

   void HomogeneousEqualityConstrainedModel::evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const {
      this->model.evaluate_objective_gradient(x, gradient);
   }

   void HomogeneousEqualityConstrainedModel::compute_jacobian_sparsity(uno_int* row_indices, uno_int* column_indices,
         uno_int row_offset, uno_int column_offset, uno_int solver_indexing, MatrixOrder matrix_order) const {
      this->model.compute_jacobian_sparsity(row_indices, column_indices, row_offset, column_offset, solver_indexing, matrix_order);

      // add the slack contributions
      size_t nonzero_index = this->model.number_jacobian_nonzeros();
      for (const auto [constraint_index, slack_index]: this->get_slacks()) {
         row_indices[nonzero_index] = static_cast<int>(constraint_index) + row_offset + solver_indexing;
         column_indices[nonzero_index] = static_cast<int>(slack_index) + column_offset + solver_indexing;
         ++nonzero_index;
      }
   }

   void HomogeneousEqualityConstrainedModel::compute_hessian_sparsity(uno_int* row_indices, uno_int* column_indices,
         uno_int solver_indexing) const {
      this->model.compute_hessian_sparsity(row_indices, column_indices, solver_indexing);
   }

   void HomogeneousEqualityConstrainedModel::evaluate_jacobian(const Vector<double>& x, double* jacobian_values) const {
      this->model.evaluate_jacobian(x, jacobian_values);

      // add the slack contributions
      size_t nonzero_index = this->model.number_jacobian_nonzeros();
      for ([[maybe_unused]] const auto _: this->get_slacks()) {
         jacobian_values[nonzero_index] = -1.;
         ++nonzero_index;
      }
   }

   void HomogeneousEqualityConstrainedModel::evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier,
         const Vector<double>& multipliers, double* hessian_values) const {
      this->model.evaluate_lagrangian_hessian(x, objective_multiplier, multipliers, hessian_values);
   }

   // linear operators for Jacobian-, Jacobian^T-, and Hessian-vector products
   void HomogeneousEqualityConstrainedModel::compute_jacobian_vector_product(const double* x, const double* vector, double* result) const {
      this->model.compute_jacobian_vector_product(x, vector, result);

      // add the slack contributions
      for (const auto [constraint_index, slack_variable_index]: this->get_slacks()) {
         result[constraint_index] -= vector[slack_variable_index];
      }
   }

   void HomogeneousEqualityConstrainedModel::compute_jacobian_transposed_vector_product(const double* x, const double* vector,
         double* result) const {
      this->model.compute_jacobian_transposed_vector_product(x, vector, result);

      // add the slack contributions
      for (const auto [constraint_index, slack_variable_index]: this->get_slacks()) {
         result[slack_variable_index] = -vector[constraint_index];
      }
   }

   void HomogeneousEqualityConstrainedModel::compute_hessian_vector_product(const double* x, const double* vector,
         double objective_multiplier, const Vector<double>& multipliers, double* result) const {
      this->model.compute_hessian_vector_product(x, vector, objective_multiplier, multipliers, result);
   }

   const std::vector<double>& HomogeneousEqualityConstrainedModel::get_variables_lower_bounds() const {
      return this->variables_lower_bounds;
   }
   
   const std::vector<double>& HomogeneousEqualityConstrainedModel::get_variables_upper_bounds() const {
      return this->variables_upper_bounds;
   }

   const SparseVector<size_t>& HomogeneousEqualityConstrainedModel::get_slacks() const {
      return this->slacks;
   }

   const Vector<size_t>& HomogeneousEqualityConstrainedModel::get_fixed_variables() const {
      return this->model.get_fixed_variables();
   }

   const std::vector<double>& HomogeneousEqualityConstrainedModel::get_constraints_lower_bounds() const {
      return this->constraints_lower_bounds;
   }

   const std::vector<double>& HomogeneousEqualityConstrainedModel::get_constraints_upper_bounds() const {
      return this->constraints_upper_bounds;
   }

   const Collection<size_t>& HomogeneousEqualityConstrainedModel::get_equality_constraints() const {
      return this->equality_constraints;
   }

   const Collection<size_t>& HomogeneousEqualityConstrainedModel::get_inequality_constraints() const {
      return this->inequality_constraints;
   }

   const Collection<size_t>& HomogeneousEqualityConstrainedModel::get_linear_constraints() const {
      return this->model.get_linear_constraints();
   }

   const Collection<size_t>& HomogeneousEqualityConstrainedModel::get_nonlinear_constraints() const {
      return this->model.get_nonlinear_constraints();
   }

   void HomogeneousEqualityConstrainedModel::initial_primal_point(Vector<double>& x) const {
      this->model.initial_primal_point(x);
      // set the slacks
      for (const auto [_, slack_index]: this->get_slacks()) {
         x[slack_index] = 0.;
      }
   }

   void HomogeneousEqualityConstrainedModel::initial_dual_point(Vector<double>& multipliers) const {
      this->model.initial_dual_point(multipliers);
   }

   void HomogeneousEqualityConstrainedModel::postprocess_solution(Iterate& iterate) const {
      // discard the slacks
      iterate.number_variables = this->model.number_variables;
      this->model.postprocess_solution(iterate);
   }

   size_t HomogeneousEqualityConstrainedModel::number_jacobian_nonzeros() const {
      return this->model.number_jacobian_nonzeros() + this->slacks.size();
   }

   size_t HomogeneousEqualityConstrainedModel::number_hessian_nonzeros() const {
      return this->model.number_hessian_nonzeros();
   }

   size_t HomogeneousEqualityConstrainedModel::number_model_objective_evaluations() const {
      return this->model.number_model_objective_evaluations();
   }

   size_t HomogeneousEqualityConstrainedModel::number_model_constraints_evaluations() const {
      return this->model.number_model_constraints_evaluations();
   }

   size_t HomogeneousEqualityConstrainedModel::number_model_objective_gradient_evaluations() const {
      return this->model.number_model_objective_gradient_evaluations();
   }

   size_t HomogeneousEqualityConstrainedModel::number_model_jacobian_evaluations() const {
      return this->model.number_model_jacobian_evaluations();
   }

   size_t HomogeneousEqualityConstrainedModel::number_model_hessian_evaluations() const {
      return this->model.number_model_hessian_evaluations();
   }

   void HomogeneousEqualityConstrainedModel::reset_number_evaluations() const {
      this->model.reset_number_evaluations();
   }
} // namespace