// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "Subproblem.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/inertia_correction_strategies/InertiaCorrectionStrategy.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/Direction.hpp"
#include "optimization/EvaluationSpace.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/ScalarMultiple.hpp"
#include "symbolic/Sum.hpp"

namespace uno {
   Subproblem::Subproblem(const OptimizationProblem& problem, Iterate& current_iterate, HessianModel& hessian_model,
      InertiaCorrectionStrategy<double>& inertia_correction_strategy):
         number_variables(problem.number_variables), number_constraints(problem.number_constraints),
         problem(problem), current_iterate(current_iterate), hessian_model(hessian_model),
         inertia_correction_strategy(inertia_correction_strategy) {
   }

   void Subproblem::compute_constraint_jacobian_sparsity(uno_int *row_indices, uno_int *column_indices, uno_int solver_indexing,
         MatrixOrder matrix_order) const {
      this->problem.compute_constraint_jacobian_sparsity(row_indices, column_indices, solver_indexing, matrix_order);
   }

   void Subproblem::compute_regularized_hessian_sparsity(uno_int *row_indices, uno_int *column_indices, uno_int solver_indexing) const {
      // sparsity of original Lagrangian Hessian
      this->problem.compute_hessian_sparsity(this->hessian_model, row_indices, column_indices, solver_indexing);

      // regularize the Hessian only if required (diagonal regularization)
      if (!this->hessian_model.is_positive_definite() && this->inertia_correction_strategy.performs_primal_regularization()) {
         size_t current_index = this->number_hessian_nonzeros();
         for (size_t variable_index: this->get_primal_regularization_variables()) {
            row_indices[current_index ] = static_cast<int>(variable_index) + solver_indexing;
            column_indices[current_index] = static_cast<int>(variable_index) + solver_indexing;
            ++current_index;
         }
      }
   }

   // lower triangular part of the symmetric augmented matrix
   void Subproblem::compute_regularized_augmented_matrix_sparsity(uno_int *row_indices, uno_int *column_indices,
         const uno_int *jacobian_row_indices, const uno_int *jacobian_column_indices, uno_int solver_indexing) const {
      // sparsity of original Lagrangian Hessian in the (1, 1) block
      this->problem.compute_hessian_sparsity(this->hessian_model, row_indices, column_indices, solver_indexing);

      // copy Jacobian of general constraints into the (2, 1) block
      const size_t hessian_offset = this->number_hessian_nonzeros();
      for (size_t nonzero_index: Range(this->number_jacobian_nonzeros())) {
         // to get the transpose, switch the order of the row/column vectors
         row_indices[hessian_offset + nonzero_index] = static_cast<int>(this->problem.number_variables) +
            jacobian_row_indices[nonzero_index] + solver_indexing;
         column_indices[hessian_offset + nonzero_index] = jacobian_column_indices[nonzero_index] + solver_indexing;
      }

      // regularize the augmented matrix only if required (diagonal regularization)
      size_t current_index = hessian_offset + this->problem.number_jacobian_nonzeros();
      if (!this->hessian_model.is_positive_definite() && this->inertia_correction_strategy.performs_primal_regularization()) {
         for (size_t variable_index: this->get_primal_regularization_variables()) {
            row_indices[current_index] = static_cast<int>(variable_index) + solver_indexing;
            column_indices[current_index] = static_cast<int>(variable_index) + solver_indexing;
            ++current_index;
         }
      }
      if (this->inertia_correction_strategy.performs_dual_regularization()) {
         for (size_t constraint_index: this->get_dual_regularization_constraints()) {
            const int shifted_constraint_index = static_cast<int>(this->number_variables + constraint_index);
            row_indices[current_index] = shifted_constraint_index + solver_indexing;
            column_indices[current_index] = shifted_constraint_index + solver_indexing;
            ++current_index;
         }
      }
   }

   void Subproblem::evaluate_lagrangian_hessian(Statistics& statistics, double* hessian_values) const {
      // evaluate the Lagrangian Hessian of the problem at the current primal-dual point
      this->problem.evaluate_lagrangian_hessian(statistics, this->hessian_model, this->current_iterate.primals,
         this->current_iterate.multipliers, hessian_values);
   }

   void Subproblem::regularize_lagrangian_hessian(Statistics& statistics, double* hessian_values) const {
      // regularize the Hessian only if necessary
      if (!this->hessian_model.is_positive_definite() && this->inertia_correction_strategy.performs_primal_regularization()) {
         const Inertia expected_inertia{this->problem.get_number_original_variables(), 0,
            this->problem.number_variables - this->problem.get_number_original_variables()};
         const size_t offset = this->number_hessian_nonzeros();
         double* primal_regularization_values = hessian_values + offset;
         this->inertia_correction_strategy.regularize_hessian(statistics, *this, hessian_values, expected_inertia,
            primal_regularization_values);
      }
   }

   void Subproblem::compute_hessian_vector_product(const double* x, const double* vector, double* result) const {
      // unregularized Hessian-vector product
      this->problem.compute_hessian_vector_product(this->hessian_model, x, vector, this->current_iterate.multipliers, result);

      // contribution of the regularization strategy
      const double regularization_factor = this->inertia_correction_strategy.get_primal_regularization_factor();
      if (0. < regularization_factor) {
         for (size_t variable_index: this->get_primal_regularization_variables()) {
            result[variable_index] += regularization_factor*vector[variable_index];
         }
      }
   }

   void Subproblem::assemble_augmented_matrix(Statistics& statistics, double* augmented_matrix_values) const {
      // evaluate the Lagrangian Hessian of the problem at the current primal-dual point
      this->problem.evaluate_lagrangian_hessian(statistics, this->hessian_model, this->current_iterate.primals,
         this->current_iterate.multipliers, augmented_matrix_values);

      // Jacobian of general constraints
      this->problem.evaluate_constraint_jacobian(this->current_iterate, augmented_matrix_values + this->number_hessian_nonzeros());
   }

   void Subproblem::regularize_augmented_matrix(Statistics& statistics, double* augmented_matrix_values,
         double dual_regularization_parameter, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver) const {
      if ((!this->hessian_model.is_positive_definite() && this->inertia_correction_strategy.performs_dual_regularization()) ||
            this->inertia_correction_strategy.performs_dual_regularization()) {
         const Inertia expected_inertia{this->number_variables, this->number_constraints, 0};

         const size_t offset = this->number_hessian_nonzeros() + this->problem.number_jacobian_nonzeros();
         double* primal_regularization_values = augmented_matrix_values + offset;
         double* dual_regularization_values = augmented_matrix_values + offset + this->get_primal_regularization_variables().size();
         this->inertia_correction_strategy.regularize_augmented_matrix(statistics, *this, augmented_matrix_values,
            dual_regularization_parameter, expected_inertia, linear_solver, primal_regularization_values, dual_regularization_values);
      }
      else {
         linear_solver.do_numerical_factorization(augmented_matrix_values);
      }
   }

   void Subproblem::assemble_primal_dual_direction(const Vector<double>& solution, Direction& direction) const {
      this->problem.assemble_primal_dual_direction(this->current_iterate, solution, direction);
   }

   void Subproblem::set_variables_bounds(std::vector<double>& variables_lower_bounds, std::vector<double>& variables_upper_bounds,
         double trust_region_radius) const {
      // bounds of original variables intersected with trust region
      for (size_t variable_index: Range(this->problem.get_number_original_variables())) {
         variables_lower_bounds[variable_index] = std::max(-trust_region_radius,
            this->problem.variable_lower_bound(variable_index) - this->current_iterate.primals[variable_index]);
         variables_upper_bounds[variable_index] = std::min(trust_region_radius,
            this->problem.variable_upper_bound(variable_index) - this->current_iterate.primals[variable_index]);
      }
      // bounds of additional variables (no trust region!)
      for (size_t variable_index: Range(this->problem.get_number_original_variables(), this->problem.number_variables)) {
         variables_lower_bounds[variable_index] = this->problem.variable_lower_bound(variable_index) - this->current_iterate.primals[variable_index];
         variables_upper_bounds[variable_index] = this->problem.variable_upper_bound(variable_index) - this->current_iterate.primals[variable_index];
      }
   }

   bool Subproblem::is_hessian_positive_definite() const {
      return this->hessian_model.is_positive_definite();
   }

   bool Subproblem::has_hessian_operator() const {
      return this->hessian_model.has_hessian_operator();
   }

   bool Subproblem::has_hessian_matrix() const {
      return this->hessian_model.has_hessian_matrix();
   }

   // two sources of curvature: the problem and the regularization strategy
   bool Subproblem::has_curvature() const {
      // the problem may have some curvature
      if (this->problem.has_curvature(this->hessian_model)) {
         return true;
      }
      else {
         // otherwise, the regularization strategy may introduce curvature
         if (!this->hessian_model.is_positive_definite() && this->inertia_correction_strategy.performs_primal_regularization()) {
            return !this->problem.get_primal_regularization_variables().empty();
         }
         return false;
      }
   }

   bool Subproblem::performs_primal_regularization() const {
      return this->inertia_correction_strategy.performs_primal_regularization();
   }

   bool Subproblem::performs_dual_regularization() const {
      return this->inertia_correction_strategy.performs_dual_regularization();
   }

   const Collection<size_t>& Subproblem::get_primal_regularization_variables() const {
      if (!this->hessian_model.is_positive_definite() && this->inertia_correction_strategy.performs_primal_regularization()) {
         return this->problem.get_primal_regularization_variables();
      }
      return this->empty_set;
   }

   const Collection<size_t>& Subproblem::get_dual_regularization_constraints() const {
      return this->problem.get_dual_regularization_constraints();
   }

   size_t Subproblem::number_jacobian_nonzeros() const {
      return this->problem.number_jacobian_nonzeros();
   }

   size_t Subproblem::number_hessian_nonzeros() const {
      return this->problem.number_hessian_nonzeros(this->hessian_model);
   }

   size_t Subproblem::number_regularized_hessian_nonzeros() const {
      size_t number_nonzeros = this->number_hessian_nonzeros();
      if (!this->hessian_model.is_positive_definite() && this->performs_primal_regularization()) {
         number_nonzeros += this->get_primal_regularization_variables().size();
      }
      return number_nonzeros;
   }

   size_t Subproblem::number_regularized_augmented_system_nonzeros() const {
      size_t number_nonzeros = this->number_hessian_nonzeros() + this->problem.number_jacobian_nonzeros();
      if (!this->hessian_model.is_positive_definite() && this->performs_primal_regularization()) {
         number_nonzeros += this->get_primal_regularization_variables().size();
      }
      if (this->performs_dual_regularization()) {
         number_nonzeros += this->get_dual_regularization_constraints().size();
      }
      return number_nonzeros;
   }

   double Subproblem::dual_regularization_factor() const {
      return this->problem.dual_regularization_factor();
   }

   // local models of progress measures
   double Subproblem::compute_predicted_infeasibility_reduction(const EvaluationSpace& evaluation_space, const Model& model,
         const Vector<double>& primal_direction, double step_length, Norm norm) const {
      // predicted infeasibility reduction: "‖c(x)‖ - ‖c(x) + ∇c(x)^T (αd)‖"
      const double current_constraint_violation = model.constraint_violation(this->current_iterate.evaluations.constraints, norm);
      Vector<double> result(model.number_constraints);
      evaluation_space.compute_constraint_jacobian_vector_product(primal_direction, result);
      const double trial_linearized_constraint_violation = model.constraint_violation(this->current_iterate.evaluations.constraints +
         step_length * result, norm);
      return current_constraint_violation - trial_linearized_constraint_violation;
   }

   std::function<double(double)> Subproblem::compute_predicted_objective_reduction(const EvaluationSpace& evaluation_space,
         const Vector<double>& primal_direction, double step_length) const {
      // predicted objective reduction: "-∇f(x)^T (αd) - α^2/2 d^T H d"
      const double directional_derivative = dot(primal_direction, this->current_iterate.evaluations.objective_gradient);
      // if the regularized Hessian is positive definite (as it usually is in line-search methods), we can compute the
      // predicted reduction with only first-order information (the directional derivative)
      const bool is_regularized_hessian_positive_definite = this->hessian_model.is_positive_definite() && this->performs_primal_regularization();
      const double quadratic_term = is_regularized_hessian_positive_definite ? 0. :
         evaluation_space.compute_hessian_quadratic_product(*this, primal_direction);
      return [=](double objective_multiplier) {
         return step_length * (-objective_multiplier*directional_derivative) - step_length*step_length/2. * quadratic_term;
      };
   }

   ProgressMeasures Subproblem::compute_predicted_reductions(const EvaluationSpace& evaluation_space, const Direction& direction,
         double step_length, Norm norm) const {
      return {
         this->compute_predicted_infeasibility_reduction(evaluation_space, this->problem.model, direction.primals,
            step_length, norm),
         this->compute_predicted_objective_reduction(evaluation_space, direction.primals, step_length),
         this->problem.compute_predicted_auxiliary_reduction(this->current_iterate, direction.primals, step_length)
      };
   }
} // namespace