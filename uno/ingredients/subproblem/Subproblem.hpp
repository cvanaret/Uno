// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SUBPROBLEM_H
#define UNO_SUBPROBLEM_H

#include "ingredients/globalization_strategies/ProgressMeasures.hpp"
#include "linear_algebra/Vector.hpp"
#include "linear_algebra/View.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "symbolic/IntegerRange.hpp"

namespace uno {
   // forward declarations
   template <typename ElementType>
   class DirectSymmetricIndefiniteLinearSolver;
   class SolverWorkspace;
   class HessianModel;
   class InertiaCorrectionStrategy;
   class Iterate;
   class Model;
   class Statistics;

   class Subproblem {
   public:
      const size_t number_variables, number_constraints;

      Subproblem(const OptimizationProblem& problem, HessianModel& hessian_model, InertiaCorrectionStrategy& inertia_correction_strategy);

      // sparsity patterns
      [[nodiscard]] View<const uno_int> get_jacobian_row_indices() const;
      [[nodiscard]] View<const uno_int> get_jacobian_column_indices() const;
      void compute_regularized_hessian_sparsity(View<uno_int> row_indices, View<uno_int> column_indices, uno_int solver_indexing) const;
      void compute_regularized_augmented_matrix_sparsity(View<uno_int> row_indices, View<uno_int> column_indices,
         uno_int solver_indexing) const;

      void evaluate_jacobian(const Iterate& current_iterate, View<double> jacobian_values, Evaluations& evaluations) const;

      // regularized Hessian
      void evaluate_lagrangian_hessian(Statistics& statistics, const Iterate& current_iterate, View<double> hessian_values) const;
      void regularize_lagrangian_hessian(Statistics& statistics, View<double> hessian_values) const;
      void compute_hessian_vector_product(const Iterate& current_iterate, View<const double> x, View<const double> vector,
         View<double> result) const;

      // augmented system
      void regularize_augmented_matrix(Statistics& statistics, View<double> primal_inertia_correction_block,
         View<double> dual_inertia_correction_block, double dual_regularization_parameter,
         DirectSymmetricIndefiniteLinearSolver<double>& linear_solver) const;
      void assemble_augmented_rhs(const Iterate& current_iterate, Evaluations& evaluations, Vector<double>& rhs) const;
      void assemble_primal_dual_direction(const Iterate& current_iterate, const Vector<double>& solution, Direction& direction) const;

      // variables bounds
      void set_variables_bounds(const Iterate& current_iterate, std::vector<double>& variables_lower_bounds, std::vector<double>& variables_upper_bounds,
         double trust_region_radius) const;

      // constraints bounds
      template <typename Array>
      void set_constraints_bounds(Array& constraints_lower_bounds, Array& constraints_upper_bounds,
         Vector<double>& constraints) const;

      [[nodiscard]] bool is_hessian_positive_definite() const;
      [[nodiscard]] bool has_hessian_operator() const;
      [[nodiscard]] bool has_hessian_matrix() const;
      [[nodiscard]] bool has_curvature() const;
      [[nodiscard]] bool has_inequality_constraints() const;
      [[nodiscard]] bool has_bound_constraints() const;

      [[nodiscard]] bool performs_primal_regularization() const;
      [[nodiscard]] bool performs_dual_regularization() const;

      [[nodiscard]] const Collection<size_t>& get_primal_regularization_variables() const;
      [[nodiscard]] const Collection<size_t>& get_dual_regularization_constraints() const;

      [[nodiscard]] size_t number_jacobian_nonzeros() const;
      [[nodiscard]] size_t number_hessian_nonzeros() const;
      [[nodiscard]] size_t number_regularized_hessian_nonzeros() const;
      [[nodiscard]] size_t number_primal_inertia_correction_nonzeros() const;
      [[nodiscard]] size_t number_dual_inertia_correction_nonzeros() const;

      [[nodiscard]] double dual_regularization_factor() const;

      // local models of progress measures
      [[nodiscard]] ProgressMeasures compute_predicted_reductions(const Iterate& current_iterate, const Direction& direction,
         double step_length, Norm norm, Evaluations& current_evaluations, const SolverWorkspace& solver_workspace) const;

      const OptimizationProblem& problem;

   protected:
      HessianModel& hessian_model;
      InertiaCorrectionStrategy& inertia_correction_strategy;
      const IntegerRange empty_set{0};
      mutable Vector<double> objective_gradient_buffer;
   };

   template <typename Array>
   void Subproblem::set_constraints_bounds(Array& constraints_lower_bounds, Array& constraints_upper_bounds,
         Vector<double>& constraints) const {
      view(constraints_lower_bounds, 0, this->number_constraints) = this->problem.get_constraints_lower_bounds() - constraints;
      view(constraints_upper_bounds, 0, this->number_constraints) = this->problem.get_constraints_upper_bounds() - constraints;
   }
} // namespace

#endif // UNO_SUBPROBLEM_H