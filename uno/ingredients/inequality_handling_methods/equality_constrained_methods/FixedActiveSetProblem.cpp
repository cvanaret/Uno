// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "FixedActiveSetProblem.hpp"
#include "optimization/Direction.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   FixedActiveSetProblem::FixedActiveSetProblem(const OptimizationProblem& problem, const ActiveSet &active_set,
      double trust_region_radius):
         OptimizationProblem(problem.model), first_reformulation(problem), active_set(active_set),
         trust_region_radius(trust_region_radius),
         variables_lower_bounds(problem.number_variables), variables_upper_bounds(problem.number_variables),
         constraints_lower_bounds(problem.number_constraints), constraints_upper_bounds(problem.number_constraints) {
      // initialize all bounds to {-INF<double>,+INF<double>}
      for (size_t i: Range(problem.number_variables)) {
         this->variables_lower_bounds[i] = -INF<double>;
         this->variables_upper_bounds[i] = INF<double>;
      }
      // set active lower bounds as equality constraints
      for (size_t variable_index: this->active_set.bounds.at_lower_bound) {
         const double lb = problem.variable_lower_bound(variable_index);
         if (-this->trust_region_radius + this->activity_tolerance <= lb) {
            this->variables_lower_bounds[variable_index] = this->variables_upper_bounds[variable_index] = lb;
         }
      }
      // set active lower bounds as equality constraints
      for (size_t variable_index: this->active_set.bounds.at_upper_bound) {
         const double ub = problem.variable_upper_bound(variable_index);
         if (ub <= this->trust_region_radius - this->activity_tolerance) {
            this->variables_lower_bounds[variable_index] = this->variables_upper_bounds[variable_index] = ub;
         }
      }

      // initialize all constraint bounds to {-INF<double>,+INF<double>}
      for (size_t constraint_index: Range(problem.number_constraints)) {
         this->constraints_lower_bounds[constraint_index] = -INF<double>;
         this->constraints_upper_bounds[constraint_index] = INF<double>;
      }
      // set active lower constraint bounds as equality constraints
      for (size_t constraint_index: this->active_set.constraints.at_lower_bound) {
         this->constraints_lower_bounds[constraint_index] = this->constraints_upper_bounds[constraint_index] =
            problem.constraint_lower_bound(constraint_index);
      }
      // set active upper constraint bounds as equality constraints
      for (size_t constraint_index: this->active_set.constraints.at_upper_bound) {
         this->constraints_lower_bounds[constraint_index] = this->constraints_upper_bounds[constraint_index] =
            problem.constraint_upper_bound(constraint_index);
      }
   }

   double FixedActiveSetProblem::variable_lower_bound(size_t variable_index) const {
      return this->variables_lower_bounds[variable_index];
   }

   double FixedActiveSetProblem::variable_upper_bound(size_t variable_index) const {
      return this->variables_upper_bounds[variable_index];
   }

   double FixedActiveSetProblem::constraint_lower_bound(size_t constraint_index) const {
      return this->constraints_lower_bounds[constraint_index];
   }

   double FixedActiveSetProblem::constraint_upper_bound(size_t constraint_index) const {
      return this->constraints_upper_bounds[constraint_index];
   }
} // namespace