// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SUBPROBLEM_H
#define UNO_SUBPROBLEM_H

#include <functional>
#include "optimization/Iterate.hpp"
#include "linear_algebra/Matrix.hpp"
#include "linear_algebra/MatrixOrder.hpp"
#include "optimization/Multipliers.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "symbolic/Range.hpp"
#include "symbolic/UnaryNegation.hpp"
#include "symbolic/VectorView.hpp"
#include "tools/Logger.hpp"

namespace uno {
   // forward declarations
   template <typename ElementType>
   class DirectSymmetricIndefiniteLinearSolver;
   class EvaluationSpace;
   class HessianModel;
   template <typename ElementType>
   class InertiaCorrectionStrategy;
   class Statistics;
   template <typename ElementType>
   class Vector;

   class Subproblem {
   public:
      const size_t number_variables, number_constraints;

      Subproblem(const OptimizationProblem& problem, Iterate& current_iterate, HessianModel& hessian_model,
         InertiaCorrectionStrategy<double>& inertia_correction_strategy);

      // sparsity patterns
      void compute_constraint_jacobian_sparsity(uno_int* row_indices, uno_int* column_indices, uno_int solver_indexing,
         MatrixOrder matrix_order) const;
      void compute_regularized_hessian_sparsity(uno_int* row_indices, uno_int* column_indices, uno_int solver_indexing) const;
      void compute_regularized_augmented_matrix_sparsity(uno_int* row_indices, uno_int* column_indices,
         const uno_int* jacobian_row_indices, const uno_int* jacobian_column_indices, uno_int solver_indexing) const;

      // regularized Hessian
      void evaluate_lagrangian_hessian(Statistics& statistics, double* hessian_values) const;
      void regularize_lagrangian_hessian(Statistics& statistics, double* hessian_values) const;
      void compute_hessian_vector_product(const double* x, const double* vector, double* result) const;

      // augmented system
      void assemble_augmented_matrix(Statistics& statistics, double* augmented_matrix_values) const;
      void regularize_augmented_matrix(Statistics& statistics, double* augmented_matrix_values,
         double dual_regularization_parameter, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver) const;
      template <typename IndexType>
      void assemble_augmented_rhs(const Vector<double>& objective_gradient, const Vector<double>& constraints,
         const Matrix<IndexType>& constraint_jacobian, Vector<double>& rhs) const;
      void assemble_primal_dual_direction(const Vector<double>& solution, Direction& direction) const;

      // variables bounds
      void set_variables_bounds(std::vector<double>& variables_lower_bounds, std::vector<double>& variables_upper_bounds,
         double trust_region_radius) const;

      // constraints bounds
      template <typename Array>
      void set_constraints_bounds(Array& constraints_lower_bounds, Array& constraints_upper_bounds,
         Vector<double>& constraints) const;

      [[nodiscard]] bool is_hessian_positive_definite() const;
      [[nodiscard]] bool has_hessian_operator() const;
      [[nodiscard]] bool has_hessian_matrix() const;
      [[nodiscard]] bool has_curvature() const;

      [[nodiscard]] bool performs_primal_regularization() const;
      [[nodiscard]] bool performs_dual_regularization() const;

      [[nodiscard]] const Collection<size_t>& get_primal_regularization_variables() const;
      [[nodiscard]] const Collection<size_t>& get_dual_regularization_constraints() const;

      [[nodiscard]] size_t number_jacobian_nonzeros() const;
      [[nodiscard]] size_t number_hessian_nonzeros() const;
      [[nodiscard]] size_t number_regularized_hessian_nonzeros() const;
      [[nodiscard]] size_t number_regularized_augmented_system_nonzeros() const;

      [[nodiscard]] double dual_regularization_factor() const;

      // local models of progress measures
      [[nodiscard]] double compute_predicted_infeasibility_reduction(const EvaluationSpace& evaluation_space, const Model& model,
         const Vector<double>& primal_direction, double step_length, Norm norm) const;
      [[nodiscard]] std::function<double(double)> compute_predicted_objective_reduction(const EvaluationSpace& evaluation_space,
         const Vector<double>& primal_direction, double step_length) const;
      [[nodiscard]] ProgressMeasures compute_predicted_reductions(const EvaluationSpace& evaluation_space,
         const Direction& direction, double step_length, Norm norm) const;

      const OptimizationProblem& problem;
      Iterate& current_iterate;

   protected:
      HessianModel& hessian_model;
      InertiaCorrectionStrategy<double>& inertia_correction_strategy;
      const ForwardRange empty_set{0};
   };

   template <typename IndexType>
   void Subproblem::assemble_augmented_rhs(const Vector<double>& objective_gradient, const Vector<double>& constraints,
         const Matrix<IndexType>& constraint_jacobian, Vector<double>& rhs) const {
      rhs.fill(0.);

      // objective gradient
      view(rhs, 0, this->number_variables) = -objective_gradient;

      // Jacobian
      for (size_t nonzero_index: Range(this->number_jacobian_nonzeros())) {
         const auto [constraint_index, variable_index, derivative] = constraint_jacobian[nonzero_index];
         rhs[static_cast<size_t>(variable_index)] +=
            this->current_iterate.multipliers.constraints[static_cast<size_t>(constraint_index)] * derivative;
      }
      // constraints
      for (size_t constraint_index: Range(this->number_constraints)) {
         rhs[this->number_variables + constraint_index] = -constraints[constraint_index];
      }
      DEBUG2 << "RHS: "; print_vector(DEBUG2, view(rhs, 0, this->number_variables + this->number_constraints)); DEBUG << '\n';
   }

   template <typename Array>
   void Subproblem::set_constraints_bounds(Array& constraints_lower_bounds, Array& constraints_upper_bounds,
         Vector<double>& constraints) const {
      for (size_t constraint_index: Range(this->problem.number_constraints)) {
         constraints_lower_bounds[constraint_index] = this->problem.constraint_lower_bound(constraint_index) - constraints[constraint_index];
         constraints_upper_bounds[constraint_index] = this->problem.constraint_upper_bound(constraint_index) - constraints[constraint_index];
      }
   }
} // namespace

#endif // UNO_SUBPROBLEM_H