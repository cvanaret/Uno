// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Subproblem.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/regularization_strategies/RegularizationStrategy.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/UnaryNegation.hpp"
#include "symbolic/VectorView.hpp"
#include "tools/Logger.hpp"

namespace uno {
   Subproblem::Subproblem(const OptimizationProblem& problem, Iterate& current_iterate, HessianModel& hessian_model,
      RegularizationStrategy<double>& regularization_strategy, double trust_region_radius):
         number_variables(problem.number_variables), number_constraints(problem.number_constraints),
         problem(problem), current_iterate(current_iterate), hessian_model(hessian_model),
         regularization_strategy(regularization_strategy), trust_region_radius(trust_region_radius) {
   }

   void Subproblem::evaluate_objective_gradient(Vector<double>& linear_objective) const {
      this->problem.evaluate_objective_gradient(this->current_iterate, linear_objective);
   }

   void Subproblem::evaluate_constraints(std::vector<double>& constraints) const {
      this->problem.evaluate_constraints(this->current_iterate, constraints);
   }

   void Subproblem::evaluate_jacobian(Vector<double>& jacobian_values) const {
      this->problem.evaluate_constraint_jacobian(this->current_iterate, jacobian_values);
   }

   void Subproblem::compute_regularized_hessian_structure(Vector<size_t>& row_indices, Vector<size_t>& column_indices) const {
      // structure of original Lagrangian Hessian
      this->problem.compute_hessian_structure(this->hessian_model, row_indices, column_indices);

      // regularize the Hessian only if required (diagonal regularization)
      if (!this->hessian_model.is_positive_definite() && this->regularization_strategy.performs_primal_regularization()) {
         for (size_t variable_index: this->get_primal_regularization_variables()) {
            row_indices.emplace_back(variable_index);
            column_indices.emplace_back(variable_index);
         }
      }
   }

   void Subproblem::compute_regularized_hessian(Statistics& statistics, Vector<double>& hessian_values) const {
      // evaluate the Lagrangian Hessian of the problem at the current primal-dual point
      this->problem.evaluate_lagrangian_hessian(statistics, this->hessian_model, this->current_iterate.primals,
         this->current_iterate.multipliers, hessian_values);

      // regularize the Hessian only if necessary
      if (!this->hessian_model.is_positive_definite() && this->regularization_strategy.performs_primal_regularization()) {
         const Inertia expected_inertia{this->problem.get_number_original_variables(), 0,
            this->problem.number_variables - this->problem.get_number_original_variables()};
         this->regularization_strategy.regularize_hessian(statistics, *this, hessian_values, expected_inertia);
      }
   }

   void Subproblem::compute_hessian_vector_product(const double* vector, double* result) const {
      // unregularized Hessian-vector product
      this->problem.compute_hessian_vector_product(this->hessian_model, vector, this->current_iterate.multipliers, result);

      // contribution of the regularization strategy
      const double regularization_factor = this->regularization_strategy.get_primal_regularization_factor();
      if (0. < regularization_factor) {
         for (size_t variable_index: this->get_primal_regularization_variables()) {
            result[variable_index] += regularization_factor*vector[variable_index];
         }
      }
   }

   void Subproblem::compute_regularized_augmented_matrix_structure(Vector<size_t>& row_indices, Vector<size_t>& column_indices) const {
      // structure of original Lagrangian Hessian
      this->problem.compute_hessian_structure(this->hessian_model, row_indices, column_indices);

      // Jacobian of general constraints: to get the transpose, switch the order of the vectors
      this->problem.compute_jacobian_structure(column_indices, row_indices);
      // shift the column indices of the Jacobian
      const size_t column_offset = this->problem.number_hessian_nonzeros(this->hessian_model);
      for (size_t nnz_index: Range(this->problem.number_jacobian_nonzeros())) {
         column_indices[column_offset + nnz_index] += this->problem.number_variables;
      }

      // regularize the augmented matrix only if required (diagonal regularization)
      if (!this->hessian_model.is_positive_definite() && this->regularization_strategy.performs_primal_regularization()) {
         for (size_t variable_index: this->get_primal_regularization_variables()) {
            row_indices.emplace_back(variable_index);
            column_indices.emplace_back(variable_index);
         }
      }
      if (this->regularization_strategy.performs_dual_regularization()) {
         for (size_t constraint_index: this->get_dual_regularization_constraints()) {
            const size_t shifted_constraint_index = this->number_variables + constraint_index;
            row_indices.emplace_back(shifted_constraint_index);
            column_indices.emplace_back(shifted_constraint_index);
         }
      }
   }

   void Subproblem::assemble_augmented_matrix(Statistics& statistics, Vector<double>& augmented_matrix_values) const {
      // evaluate the Lagrangian Hessian of the problem at the current primal-dual point
      this->problem.evaluate_lagrangian_hessian(statistics, this->hessian_model, this->current_iterate.primals,
         this->current_iterate.multipliers, augmented_matrix_values);

      // Jacobian of general constraints
      this->problem.evaluate_constraint_jacobian(this->current_iterate, augmented_matrix_values);
   }

   void Subproblem::regularize_augmented_matrix(Statistics& statistics, Vector<double>& augmented_matrix_values,
         double dual_regularization_parameter, DirectSymmetricIndefiniteLinearSolver<size_t, double>& linear_solver) const {
      const Inertia expected_inertia{this->number_variables, this->number_constraints, 0};
      this->regularization_strategy.regularize_augmented_matrix(statistics, *this, augmented_matrix_values,
         dual_regularization_parameter, expected_inertia, linear_solver);
   }

   void Subproblem::assemble_augmented_rhs(const Vector<double>& objective_gradient, const std::vector<double>& constraints,
         RectangularMatrix<double>& constraint_jacobian, Vector<double>& rhs) const {
      rhs.fill(0.);

      // objective gradient
      view(rhs, 0, this->number_variables) = -objective_gradient;

      // constraint: evaluations and gradients
      for (size_t constraint_index: Range(this->number_constraints)) {
         // Lagrangian
         /*
         if (this->current_iterate.multipliers.constraints[constraint_index] != 0.) {
            for (const auto [variable_index, derivative]: constraint_jacobian[constraint_index]) {
               rhs[variable_index] += this->current_iterate.multipliers.constraints[constraint_index] * derivative;
            }
         }
         */
         // constraints
         rhs[this->number_variables + constraint_index] = -constraints[constraint_index];
      }
      DEBUG2 << "RHS: "; print_vector(DEBUG2, view(rhs, 0, this->number_variables + this->number_constraints)); DEBUG << '\n';
   }

   void Subproblem::assemble_primal_dual_direction(const Vector<double>& solution, Direction& direction) const {
      this->problem.assemble_primal_dual_direction(this->current_iterate, solution, direction);
   }

   void Subproblem::set_variables_bounds(std::vector<double>& variables_lower_bounds, std::vector<double>& variables_upper_bounds) const {
      // bounds of original variables intersected with trust region
      for (size_t variable_index: Range(this->problem.get_number_original_variables())) {
         variables_lower_bounds[variable_index] = std::max(-this->trust_region_radius,
            this->problem.variable_lower_bound(variable_index) - this->current_iterate.primals[variable_index]);
         variables_upper_bounds[variable_index] = std::min(this->trust_region_radius,
            this->problem.variable_upper_bound(variable_index) - this->current_iterate.primals[variable_index]);
      }
      // bounds of additional variables (no trust region!)
      for (size_t variable_index: Range(this->problem.get_number_original_variables(), this->problem.number_variables)) {
         variables_lower_bounds[variable_index] = this->problem.variable_lower_bound(variable_index) - this->current_iterate.primals[variable_index];
         variables_upper_bounds[variable_index] = this->problem.variable_upper_bound(variable_index) - this->current_iterate.primals[variable_index];
      }
   }

   bool Subproblem::has_implicit_hessian_representation() const {
      return this->hessian_model.has_implicit_representation();
   }

   bool Subproblem::has_explicit_hessian_representation() const {
      return this->hessian_model.has_explicit_representation();
   }

   // two sources of curvature: the problem and the regularization strategy
   bool Subproblem::has_curvature() const {
      // the problem may have some curvature
      if (this->problem.has_curvature(this->hessian_model)) {
         return true;
      }
      else {
         // otherwise, the regularization strategy may introduce curvature
         if (!this->hessian_model.is_positive_definite() && this->regularization_strategy.performs_primal_regularization()) {
            return !this->problem.get_primal_regularization_variables().empty();
         }
         return false;
      }
   }

   bool Subproblem::performs_primal_regularization() const {
      return this->regularization_strategy.performs_primal_regularization();
   }

   bool Subproblem::performs_dual_regularization() const {
      return this->regularization_strategy.performs_dual_regularization();
   }

   const Collection<size_t>& Subproblem::get_primal_regularization_variables() const {
      return this->problem.get_primal_regularization_variables();
   }

   const Collection<size_t>& Subproblem::get_dual_regularization_constraints() const {
      return this->problem.get_dual_regularization_constraints();
   }

   size_t Subproblem::number_augmented_system_nonzeros() const {
      return this->problem.number_hessian_nonzeros(this->hessian_model) + this->problem.number_jacobian_nonzeros();
   }

   size_t Subproblem::regularization_size() const {
      const size_t primal_regularization_size = this->get_primal_regularization_variables().size();
      const size_t dual_regularization_size = this->get_dual_regularization_constraints().size();
      const size_t regularization_size =
         (this->performs_primal_regularization() ? primal_regularization_size : 0) +
         (this->performs_dual_regularization() ? dual_regularization_size : 0);
      return regularization_size;
   }

   double Subproblem::dual_regularization_factor() const {
      return this->problem.dual_regularization_factor();
   }
} // namespace