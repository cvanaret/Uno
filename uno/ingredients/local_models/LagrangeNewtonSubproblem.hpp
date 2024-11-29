// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <vector>
#include "optimization/Iterate.hpp"
#include "reformulation/OptimizationProblem.hpp"
#include "symbolic/VectorView.hpp"

namespace uno {
   // forward declaration
   class HessianModel;

   class LagrangeNewtonSubproblem {
   public:
      LagrangeNewtonSubproblem(const OptimizationProblem& problem, Iterate& current_iterate, const Vector<double>& current_multipliers,
            const HessianModel& hessian_model, double trust_region_radius);

      const size_t number_variables;
      const size_t number_constraints;

      template <typename Array>
      void set_direction_bounds(Array& lower_bounds, Array& upper_bounds) const;

      template <typename Array>
      void set_constraint_bounds(const Vector<double>& current_constraints, Array& lower_bounds, Array& upper_bounds) const;

      void evaluate_objective_gradient(SparseVector<double>& gradient) const;
      void evaluate_constraints(Vector<double>& constraints) const;
      void evaluate_constraint_jacobian(RectangularMatrix<double>& jacobian) const;
      void evaluate_lagrangian_hessian(const Vector<double>& x, SymmetricMatrix<size_t, double>& hessian) const;
      void compute_hessian_vector_product(const Vector<double>& x, Vector<double>& result) const;

   protected:
      const OptimizationProblem& problem;
      Iterate& current_iterate;
      const Vector<double>& current_multipliers;
      const HessianModel& hessian_model;
      const double trust_region_radius;
   };

   template <typename Array>
   void LagrangeNewtonSubproblem::set_direction_bounds(Array& lower_bounds, Array& upper_bounds) const {
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
   void LagrangeNewtonSubproblem::set_constraint_bounds(const Vector<double>& current_constraints, Array& lower_bounds, Array& upper_bounds) const {
      for (size_t constraint_index: Range(this->problem.number_constraints)) {
         lower_bounds[constraint_index] = this->problem.constraint_lower_bound(constraint_index) - current_constraints[constraint_index];
         upper_bounds[constraint_index] = this->problem.constraint_upper_bound(constraint_index) - current_constraints[constraint_index];
      }
   }
} // namespace