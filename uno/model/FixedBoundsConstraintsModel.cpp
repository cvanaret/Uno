// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "FixedBoundsConstraintsModel.hpp"
#include "linear_algebra/VectorView.hpp"
#include "optimization/Iterate.hpp"

namespace uno {
   FixedBoundsConstraintsModel::FixedBoundsConstraintsModel(const Model& original_model) :
         Model(original_model.name + " -> no fixed bounds", original_model.number_variables,
            // move the fixed variables to the set of general constraints
            original_model.number_constraints + original_model.get_fixed_variables().size(),
            original_model.optimization_sense, original_model.lagrangian_sign_convention),
         model(original_model),
         equality_constraints(concatenate(this->model.get_equality_constraints(), Range(this->model.number_constraints, this->number_constraints))),
         linear_constraints(concatenate(this->model.get_linear_constraints(), Range(this->model.number_constraints, this->number_constraints))),
         variables_lower_bounds(original_model.get_variables_lower_bounds()),
         variables_upper_bounds(original_model.get_variables_upper_bounds()),
         constraints_lower_bounds(this->number_constraints),
         constraints_upper_bounds(this->number_constraints) {
      // copy the original constraints bounds
      view(constraints_lower_bounds, 0, original_model.number_constraints) = original_model.get_constraints_lower_bounds();
      view(constraints_upper_bounds, 0, original_model.number_constraints) = original_model.get_constraints_upper_bounds();
      // handle the fixed variables
      size_t fixed_variable_constraint_index = original_model.number_constraints;
      for (size_t fixed_variable_index: original_model.get_fixed_variables()) {
         // relax the bounds of the fixed variables
         this->variables_lower_bounds[fixed_variable_index] = -INF<double>;
         this->variables_upper_bounds[fixed_variable_index] = INF<double>;
         const double fixed_value = original_model.get_variables_lower_bounds()[fixed_variable_index];
         // set the bounds of the corresponding new constraint
         this->constraints_lower_bounds[fixed_variable_constraint_index] = fixed_value;
         this->constraints_upper_bounds[fixed_variable_constraint_index] = fixed_value;
         ++fixed_variable_constraint_index;
      }
   }

   ProblemType FixedBoundsConstraintsModel::get_problem_type() const {
      return this->model.get_problem_type();
   }

   bool FixedBoundsConstraintsModel::has_jacobian_operator() const {
      return this->model.has_jacobian_operator();
   }

   bool FixedBoundsConstraintsModel::has_jacobian_transposed_operator() const {
      return this->model.has_jacobian_transposed_operator();
   }

   bool FixedBoundsConstraintsModel::has_hessian_operator() const {
      return this->model.has_hessian_operator();
   }

   bool FixedBoundsConstraintsModel::has_hessian_matrix() const {
      return this->model.has_hessian_matrix();
   }

   double FixedBoundsConstraintsModel::evaluate_objective(const Vector<double>& x) const {
      return this->model.evaluate_objective(x);
   }

   void FixedBoundsConstraintsModel::evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const {
      this->model.evaluate_constraints(x, constraints);
      // add the fixed variables
      size_t current_constraint = this->model.number_constraints;
      for (size_t fixed_variable_index: this->model.get_fixed_variables()) {
         constraints[current_constraint] = x[fixed_variable_index];
         ++current_constraint;
      }
   }

   void FixedBoundsConstraintsModel::evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const {
      this->model.evaluate_objective_gradient(x, gradient);
   }

   void FixedBoundsConstraintsModel::compute_jacobian_sparsity(uno_int* row_indices, uno_int* column_indices,
         uno_int solver_indexing, MatrixOrder matrix_order) const {
      // original constraints
      this->model.compute_jacobian_sparsity(row_indices, column_indices, solver_indexing, matrix_order);

      // fixed variables (as linear constraints)
      int constraint_index = static_cast<int>(this->model.number_constraints);
      size_t current_index = this->model.number_jacobian_nonzeros();
      for (size_t fixed_variable_index: this->model.get_fixed_variables()) {
         row_indices[current_index] = constraint_index + solver_indexing;
         column_indices[current_index] = static_cast<int>(fixed_variable_index) + solver_indexing;
         ++constraint_index;
         ++current_index;
      }
   }

   void FixedBoundsConstraintsModel::compute_hessian_sparsity(uno_int* row_indices, uno_int* column_indices, uno_int solver_indexing) const {
      this->model.compute_hessian_sparsity(row_indices, column_indices, solver_indexing);
   }

   void FixedBoundsConstraintsModel::evaluate_jacobian(const Vector<double>& x, double* jacobian_values) const {
      this->model.evaluate_jacobian(x, jacobian_values);

      // add the contributions of the fixed variables
      size_t nonzero_index = this->model.number_jacobian_nonzeros();
      for ([[maybe_unused]] size_t _: this->model.get_fixed_variables()) {
         jacobian_values[nonzero_index] = 1.;
         ++nonzero_index;
      }
   }

   void FixedBoundsConstraintsModel::evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier,
         const Vector<double>& multipliers, double* hessian_values) const {
      this->model.evaluate_lagrangian_hessian(x, objective_multiplier, multipliers, hessian_values);
   }

   // linear operators for Jacobian-, Jacobian^T-, and Hessian-vector products
   void FixedBoundsConstraintsModel::compute_jacobian_vector_product(const double* x, const double* vector, double* result) const {
      this->model.compute_jacobian_vector_product(x, vector, result);

      // add the contributions of the fixed variables
      size_t constraint_index = this->number_constraints;
      for (size_t fixed_variable_index: this->model.get_fixed_variables()) {
         result[constraint_index] = vector[fixed_variable_index];
         ++constraint_index;
      }
   }

   void FixedBoundsConstraintsModel::compute_jacobian_transposed_vector_product(const double* x, const double* vector,
         double* result) const {
      this->model.compute_jacobian_transposed_vector_product(x, vector, result);

      // add the contributions of the fixed variables
      size_t constraint_index = this->number_constraints;
      for (size_t fixed_variable_index: this->model.get_fixed_variables()) {
         result[fixed_variable_index] += vector[constraint_index];
         ++constraint_index;
      }
   }

   void FixedBoundsConstraintsModel::compute_hessian_vector_product(const double* x, const double* vector,
         double objective_multiplier, const Vector<double>& multipliers, double* result) const {
      this->model.compute_hessian_vector_product(x, vector, objective_multiplier, multipliers, result);
   }

   const std::vector<double>& FixedBoundsConstraintsModel::get_variables_lower_bounds() const {
      return this->variables_lower_bounds;
   }

   const std::vector<double>& FixedBoundsConstraintsModel::get_variables_upper_bounds() const {
      return this->variables_upper_bounds;
   }

   const SparseVector<size_t>& FixedBoundsConstraintsModel::get_slacks() const {
      return this->model.get_slacks();
   }

   const Vector<size_t>& FixedBoundsConstraintsModel::get_fixed_variables() const {
      return this->fixed_variables;
   }

   const std::vector<double>& FixedBoundsConstraintsModel::get_constraints_lower_bounds() const {
      return this->constraints_lower_bounds;
   }

   const std::vector<double>& FixedBoundsConstraintsModel::get_constraints_upper_bounds() const {
      return this->constraints_upper_bounds;
   }

   const Collection<size_t>& FixedBoundsConstraintsModel::get_equality_constraints() const {
      return this->equality_constraints;
   }

   const Collection<size_t>& FixedBoundsConstraintsModel::get_inequality_constraints() const {
      return this->model.get_inequality_constraints();
   }

   const Collection<size_t>& FixedBoundsConstraintsModel::get_linear_constraints() const {
      return this->linear_constraints;
   }

   const Collection<size_t>& FixedBoundsConstraintsModel::get_nonlinear_constraints() const {
      return this->model.get_nonlinear_constraints();
   }

   void FixedBoundsConstraintsModel::initial_primal_point(Vector<double>& x) const {
      this->model.initial_primal_point(x);
      // set the fixed variables
      for (size_t variable_index: this->model.get_fixed_variables()) {
         x[variable_index] = this->model.get_variables_lower_bounds()[variable_index];
      }
   }

   void FixedBoundsConstraintsModel::initial_dual_point(Vector<double>& multipliers) const {
      this->model.initial_dual_point(multipliers);
   }

   void FixedBoundsConstraintsModel::postprocess_solution(Iterate& iterate) const {
      // move the multipliers back from the general constraints to the bound constraints
      size_t current_constraint = this->model.number_constraints;
      for (size_t variable_index: this->model.get_fixed_variables()) {
         const double constraint_multiplier = iterate.multipliers.constraints[current_constraint];
         if (0. < constraint_multiplier) {
            iterate.multipliers.lower_bounds[variable_index] = constraint_multiplier;
         }
         else {
            iterate.multipliers.upper_bounds[variable_index] = constraint_multiplier;
         }
         ++current_constraint;
      }
      this->model.postprocess_solution(iterate);
   }

   size_t FixedBoundsConstraintsModel::number_jacobian_nonzeros() const {
      return this->model.number_jacobian_nonzeros() + this->model.get_fixed_variables().size();
   }

   size_t FixedBoundsConstraintsModel::number_hessian_nonzeros() const {
      return this->model.number_hessian_nonzeros();
   }

   size_t FixedBoundsConstraintsModel::number_model_objective_evaluations() const {
      return this->model.number_model_objective_evaluations();
   }

   size_t FixedBoundsConstraintsModel::number_model_constraints_evaluations() const {
      return this->model.number_model_constraints_evaluations();
   }

   size_t FixedBoundsConstraintsModel::number_model_objective_gradient_evaluations() const {
      return this->model.number_model_objective_gradient_evaluations();
   }

   size_t FixedBoundsConstraintsModel::number_model_jacobian_evaluations() const {
      return this->model.number_model_jacobian_evaluations();
   }

   size_t FixedBoundsConstraintsModel::number_model_hessian_evaluations() const {
      return this->model.number_model_hessian_evaluations();
   }

   void FixedBoundsConstraintsModel::reset_number_evaluations() const {
      this->model.reset_number_evaluations();
   }
} // namespace