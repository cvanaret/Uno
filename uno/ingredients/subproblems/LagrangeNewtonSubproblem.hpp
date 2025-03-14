// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LAGRANGENEWTONSUBPROBLEM_H
#define UNO_LAGRANGENEWTONSUBPROBLEM_H

#include <cstddef>
#include <vector>
#include "optimization/OptimizationProblem.hpp"
#include "optimization/Iterate.hpp"

namespace uno {
   // forward declarations
   template <typename IndexType, typename ElementType>
   class DirectSymmetricIndefiniteLinearSolver;
   class HessianModel;
   class Iterate;
   class Multipliers;
   class OptimizationProblem;
   template <typename ElementType>
   class RectangularMatrix;
    template <typename ElementType>
    class RegularizationStrategy;
   template <typename ElementType>
   class SparseVector;
   class Statistics;
   template <typename IndexType, typename ElementType>
   class SymmetricMatrix;
   template <typename ElementType>
   class Vector;
   class WarmstartInformation;

   class LagrangeNewtonSubproblem {
   public:
      const size_t number_variables;
      const size_t number_constraints;

      LagrangeNewtonSubproblem(const OptimizationProblem& problem, Iterate& current_iterate, const Multipliers& current_multipliers,
            HessianModel& hessian_model, RegularizationStrategy<double>& regularization_strategy, double trust_region_radius);

      void evaluate_objective_gradient(SparseVector<double>& objective_gradient);
      void evaluate_constraints(Vector<double>& constraints);
      void evaluate_constraint_jacobian(RectangularMatrix<double>& jacobian);
      void compute_lagrangian_gradient(SparseVector<double>& objective_gradient, RectangularMatrix<double>& jacobian, Vector<double>& gradient) const;
      void evaluate_hessian(SymmetricMatrix<size_t, double>& hessian);
      void evaluate_functions(SparseVector<double>& objective_gradient, Vector<double>& constraints, RectangularMatrix<double>& jacobian,
         SymmetricMatrix<size_t, double>& hessian, WarmstartInformation& warmstart_information);

      void assemble_augmented_rhs(SparseVector<double>& objective_gradient, Vector<double>& constraints, RectangularMatrix<double>& jacobian,
         Vector<double>& rhs) const;

      void finalize_augmented_matrix(Statistics& statistics, SymmetricMatrix<size_t, double>& augmented_matrix,
         DirectSymmetricIndefiniteLinearSolver<size_t, double>& linear_solver, WarmstartInformation& warmstart_information);

      template <typename Array>
      void set_variable_bounds(Array& lower_bounds, Array& upper_bounds);
      template <typename Array>
      void set_constraint_bounds(const Vector<double>& constraints, Array& lower_bounds, Array& upper_bounds);

   protected:
      const OptimizationProblem& problem;
      Iterate& current_iterate;
      const Multipliers& current_multipliers;
      HessianModel& hessian_model;
      RegularizationStrategy<double>& regularization_strategy;
      const double trust_region_radius;
   };

   template <typename Array>
   void LagrangeNewtonSubproblem::set_variable_bounds(Array& lower_bounds, Array& upper_bounds) {
      // bounds of original variables intersected with trust region
      for (size_t variable_index: Range(this->problem.get_number_original_variables())) {
         lower_bounds[variable_index] = std::max(-this->trust_region_radius,
               this->problem.variable_lower_bound(variable_index) - this->current_iterate.primals[variable_index]);
         upper_bounds[variable_index] = std::min(this->trust_region_radius,
               this->problem.variable_upper_bound(variable_index) - this->current_iterate.primals[variable_index]);
      }
      // bounds of additional variables (no trust region!)
      for (size_t variable_index: Range(this->problem.get_number_original_variables(), this->problem.number_variables)) {
         lower_bounds[variable_index] = this->problem.variable_lower_bound(variable_index) - this->current_iterate.primals[variable_index];
         upper_bounds[variable_index] = this->problem.variable_upper_bound(variable_index) - this->current_iterate.primals[variable_index];
      }
   }

   template <typename Array>
   void LagrangeNewtonSubproblem::set_constraint_bounds(const Vector<double>& constraints, Array& lower_bounds, Array& upper_bounds) {
      for (size_t constraint_index: Range(this->problem.number_constraints)) {
         lower_bounds[constraint_index] = problem.constraint_lower_bound(constraint_index) - constraints[constraint_index];
         upper_bounds[constraint_index] = problem.constraint_upper_bound(constraint_index) - constraints[constraint_index];
      }
   }
} // namespace

#endif // UNO_LAGRANGENEWTONSUBPROBLEM_H