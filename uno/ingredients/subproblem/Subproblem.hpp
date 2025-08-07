// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SUBPROBLEM_H
#define UNO_SUBPROBLEM_H

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
   template <typename IndexType, typename ElementType>
   class DirectSymmetricIndefiniteLinearSolver;
   class HessianModel;
   template <typename ElementType>
   class RegularizationStrategy;
   class Statistics;
   template <typename ElementType>
   class Vector;

   class Subproblem {
   public:
      const size_t number_variables, number_constraints;

      Subproblem(const OptimizationProblem& problem, Iterate& current_iterate, HessianModel& hessian_model,
         RegularizationStrategy<double>& regularization_strategy, double trust_region_radius);

      // sparsity patterns
      void compute_constraint_jacobian_sparsity(size_t* row_indices, size_t* column_indices, size_t solver_indexing,
         MatrixOrder matrix_order) const;
      void compute_regularized_hessian_sparsity(size_t* row_indices, size_t* column_indices, size_t solver_indexing) const;
      void compute_regularized_augmented_matrix_sparsity(size_t* row_indices, size_t* column_indices, const size_t* jacobian_row_indices,
         const size_t* jacobian_column_indices, size_t solver_indexing) const;

      // constraints, objective gradient and Jacobian
      void evaluate_objective_gradient(Vector<double>& linear_objective) const;
      void evaluate_constraints(std::vector<double>& constraints) const;
      void evaluate_constraint_jacobian(double* jacobian_values) const;

      // regularized Hessian
      void compute_regularized_hessian(Statistics& statistics, Vector<double>& /*hessian_values*/) const;
      void compute_hessian_vector_product(const double* vector, double* result) const;

      // augmented system
      void assemble_augmented_matrix(Statistics& statistics, Vector<double>& augmented_matrix_values) const;
      void regularize_augmented_matrix(Statistics& statistics, Vector<double>& augmented_matrix_values,
         double dual_regularization_parameter, DirectSymmetricIndefiniteLinearSolver<size_t, double>& linear_solver) const;
      template <typename IndexType>
      void assemble_augmented_rhs(const Vector<double>& objective_gradient, const std::vector<double>& constraints,
         const Matrix<IndexType>& constraint_jacobian, Vector<double>& rhs) const;
      void assemble_primal_dual_direction(const Vector<double>& solution, Direction& direction) const;

      // variables bounds
      void set_variables_bounds(std::vector<double>& variables_lower_bounds, std::vector<double>& variables_upper_bounds) const;

      // constraints bounds
      template <typename Array>
      void set_constraints_bounds(Array& constraints_lower_bounds, Array& constraints_upper_bounds,
         std::vector<double>& constraints) const;

      [[nodiscard]] bool has_implicit_hessian_representation() const;
      [[nodiscard]] bool has_explicit_hessian_representation() const;
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

   protected:
      const OptimizationProblem& problem;
      Iterate& current_iterate;
      HessianModel& hessian_model;
      RegularizationStrategy<double>& regularization_strategy;
      const double trust_region_radius;

      [[nodiscard]] size_t regularization_size() const;
   };

   template <typename IndexType>
   void Subproblem::assemble_augmented_rhs(const Vector<double>& objective_gradient, const std::vector<double>& constraints,
         const Matrix<IndexType>& constraint_jacobian, Vector<double>& rhs) const {
      rhs.fill(0.);

      // objective gradient
      view(rhs, 0, this->number_variables) = -objective_gradient;

      // Jacobian
      for (size_t nonzero_index: Range(this->number_jacobian_nonzeros())) {
         const auto [constraint_index, variable_index, derivative] = constraint_jacobian[nonzero_index];
         rhs[variable_index] += this->current_iterate.multipliers.constraints[constraint_index] * derivative;
      }
      // constraints
      for (size_t constraint_index: Range(this->number_constraints)) {
         rhs[this->number_variables + constraint_index] = -constraints[constraint_index];
      }
      DEBUG2 << "RHS: "; print_vector(DEBUG2, view(rhs, 0, this->number_variables + this->number_constraints)); DEBUG << '\n';
   }

   template <typename Array>
   void Subproblem::set_constraints_bounds(Array& constraints_lower_bounds, Array& constraints_upper_bounds,
         std::vector<double>& constraints) const {
      for (size_t constraint_index: Range(this->problem.number_constraints)) {
         constraints_lower_bounds[constraint_index] = this->problem.constraint_lower_bound(constraint_index) - constraints[constraint_index];
         constraints_upper_bounds[constraint_index] = this->problem.constraint_upper_bound(constraint_index) - constraints[constraint_index];
      }
   }
} // namespace

#endif // UNO_SUBPROBLEM_H