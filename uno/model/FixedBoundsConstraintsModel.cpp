// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "FixedBoundsConstraintsModel.hpp"
#include "optimization/Iterate.hpp"

namespace uno {
   FixedBoundsConstraintsModel::FixedBoundsConstraintsModel(const Model& original_model) :
         Model(original_model.name + " -> no fixed bounds", original_model.number_variables,
            // move the fixed variables to the set of general constraints
            original_model.number_constraints + original_model.get_fixed_variables().size(),
            original_model.optimization_sense),
         model(original_model),
         equality_constraints(concatenate(this->model.get_equality_constraints(), Range(this->model.number_constraints, this->number_constraints))),
         linear_constraints(concatenate(this->model.get_linear_constraints(), Range(this->model.number_constraints, this->number_constraints))) {
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

   void FixedBoundsConstraintsModel::compute_constraint_jacobian_sparsity(int* row_indices, int* column_indices,
         int solver_indexing, MatrixOrder matrix_order) const {
      // original constraints
      this->model.compute_constraint_jacobian_sparsity(row_indices, column_indices, solver_indexing, matrix_order);

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

   void FixedBoundsConstraintsModel::compute_hessian_sparsity(int* row_indices, int* column_indices, int solver_indexing) const {
      this->model.compute_hessian_sparsity(row_indices, column_indices, solver_indexing);
   }

   void FixedBoundsConstraintsModel::evaluate_constraint_jacobian(const Vector<double>& x, double* jacobian_values) const {
      this->model.evaluate_constraint_jacobian(x, jacobian_values);

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

   double FixedBoundsConstraintsModel::variable_lower_bound(size_t variable_index) const {
      if (this->model.variable_lower_bound(variable_index) == this->model.variable_upper_bound(variable_index)) {
      // remove bounds of fixed variables
         return -INF<double>;
      }
      return this->model.variable_lower_bound(variable_index);
   }

   double FixedBoundsConstraintsModel::variable_upper_bound(size_t variable_index) const {
      if (this->model.variable_lower_bound(variable_index) == this->model.variable_upper_bound(variable_index)) {
      // remove bounds of fixed variables
         return INF<double>;
      }
      return this->model.variable_upper_bound(variable_index);
   }

   const SparseVector<size_t>& FixedBoundsConstraintsModel::get_slacks() const {
      return this->model.get_slacks();
   }

   const Vector<size_t>& FixedBoundsConstraintsModel::get_fixed_variables() const {
      return this->fixed_variables;
   }

   double FixedBoundsConstraintsModel::constraint_lower_bound(size_t constraint_index) const {
      if (constraint_index < this->model.number_constraints) {
// original constraint
         return this->model.constraint_lower_bound(constraint_index);
      }
      else {
// fixed variable
         const size_t variable_index = this->model.get_fixed_variables()[constraint_index - this->model.number_constraints];
         return this->model.variable_lower_bound(variable_index);
      }
   }

   double FixedBoundsConstraintsModel::constraint_upper_bound(size_t constraint_index) const {
      if (constraint_index < this->model.number_constraints) {
         // original constraint
         return this->model.constraint_upper_bound(constraint_index);
      }
      else {
         // fixed variable
         const size_t variable_index = this->model.get_fixed_variables()[constraint_index - this->model.number_constraints];
         return this->model.variable_lower_bound(variable_index);
      }
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

   void FixedBoundsConstraintsModel::initial_primal_point(Vector<double>& x) const {
      this->model.initial_primal_point(x);
      // set the fixed variables
      for (size_t variable_index: this->model.get_fixed_variables()) {
         x[variable_index] = this->model.variable_lower_bound(variable_index);
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
} // namespace