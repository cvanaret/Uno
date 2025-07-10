// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SUBPROBLEM_H
#define UNO_SUBPROBLEM_H

#include "ingredients/regularization_strategies/RegularizationStrategy.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "symbolic/Range.hpp"

namespace uno {
   // forward declarations
   template <typename IndexType, typename ElementType>
   class DirectSymmetricIndefiniteLinearSolver;
   class HessianModel;
   class Iterate;
   class Multipliers;
   template <typename ElementType>
   class RectangularMatrix;
   template <typename ElementType>
   class SparseVector;
   class Statistics;
   template <typename IndexType, typename ElementType>
   class SymmetricMatrix;
   template <typename ElementType>
   class Vector;

   class Subproblem {
   public:
      const size_t number_variables, number_constraints;

      Subproblem(const OptimizationProblem& problem, Iterate& current_iterate, const Multipliers& current_multipliers,
         HessianModel& hessian_model, RegularizationStrategy<double>& regularization_strategy, double trust_region_radius);

      // constraints, objective gradient and Jacobian
      void evaluate_objective_gradient(SparseVector<double>& linear_objective) const;
      void evaluate_constraints(std::vector<double>& constraints) const;
      void evaluate_jacobian(RectangularMatrix<double>& constraint_jacobian) const;

      // regularized Hessian
      void compute_regularized_hessian(Statistics& statistics, SymmetricMatrix<size_t, double>& hessian) const;
      template <typename Vector1, typename Vector2>
      void compute_hessian_vector_product(const Vector1& vector, Vector2& result) const;

      // augmented system
      void assemble_augmented_matrix(Statistics& statistics, SymmetricMatrix<size_t, double>& augmented_matrix,
         RectangularMatrix<double>& constraint_jacobian) const;
      void regularize_augmented_matrix(Statistics& statistics, SymmetricMatrix<size_t, double>& augmented_matrix,
         double dual_regularization_parameter, DirectSymmetricIndefiniteLinearSolver<size_t, double>& linear_solver) const;
      void assemble_augmented_rhs(const SparseVector<double>& objective_gradient, const std::vector<double>& constraints,
         RectangularMatrix<double>& constraint_jacobian, Vector<double>& rhs) const;
      void assemble_primal_dual_direction(const Vector<double>& solution, Direction& direction) const;

      // variables bounds
      void set_variables_bounds(std::vector<double>& variables_lower_bounds, std::vector<double>& variables_upper_bounds) const;

      // constraints bounds
      template <typename Array>
      void set_constraints_bounds(Array& constraints_lower_bounds, Array& constraints_upper_bounds,
         std::vector<double>& constraints) const;

      [[nodiscard]] bool has_implicit_hessian_representation() const;
      [[nodiscard]] bool has_explicit_hessian_representation() const;

      [[nodiscard]] double dual_regularization_factor() const;

   protected:
      const OptimizationProblem& problem;
      Iterate& current_iterate;
      const Multipliers& current_multipliers;
      HessianModel& hessian_model;
      RegularizationStrategy<double>& regularization_strategy;
      const double trust_region_radius;
   };

   template <typename Array>
   void Subproblem::set_constraints_bounds(Array& constraints_lower_bounds, Array& constraints_upper_bounds,
         std::vector<double>& constraints) const {
      for (size_t constraint_index: Range(this->problem.number_constraints)) {
         constraints_lower_bounds[constraint_index] = this->problem.constraint_lower_bound(constraint_index) - constraints[constraint_index];
         constraints_upper_bounds[constraint_index] = this->problem.constraint_upper_bound(constraint_index) - constraints[constraint_index];
      }
   }

   template <typename Vector1, typename Vector2>
   void Subproblem::compute_hessian_vector_product(const Vector1& vector, Vector2& result) const {
      // copy to and from uno::Vector
      Vector<double> uno_vector(this->number_variables); // TODO wrap
      for (size_t index: Range(this->number_variables)) {
         uno_vector[index] = vector[index];
      }
      Vector<double> uno_result(this->number_variables); // TODO wrap
      // unregularized Hessian-vector product
      this->problem.compute_hessian_vector_product(this->hessian_model, uno_vector, this->current_multipliers, uno_result);
      // contribution of the regularization strategy
      const double regularization_factor = this->regularization_strategy.get_primal_regularization_factor();
      if (0. < regularization_factor) {
         for (size_t index: Range(this->number_variables)) {
            result[index] += regularization_factor*vector[index];
         }
      }
      for (size_t index: Range(this->number_variables)) {
         result[index] = uno_result[index];
      }
   }
} // namespace

#endif // UNO_SUBPROBLEM_H