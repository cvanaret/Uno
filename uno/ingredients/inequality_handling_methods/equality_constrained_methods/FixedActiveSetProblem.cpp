// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "FixedActiveSetProblem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Infinity.hpp"
#include "tools/Logger.hpp"

namespace uno {
   FixedActiveSetProblem::FixedActiveSetProblem(const OptimizationProblem& problem, const Vector<double>& direction_primals,
      const ActiveSet &active_set, double trust_region_radius):
         OptimizationProblem(problem.model), problem(problem), direction_primals(direction_primals),
         active_set(active_set),
         trust_region_radius(trust_region_radius),
         variables_lower_bounds(problem.number_variables), variables_upper_bounds(problem.number_variables),
         constraints_lower_bounds(problem.number_constraints), constraints_upper_bounds(problem.number_constraints) {
      // variables: initialize all bounds to {-inf, +inf}
      for (size_t i: Range(problem.number_variables)) {
         this->variables_lower_bounds[i] = -INF<double>;
         this->variables_upper_bounds[i] = INF<double>;
      }
      // variables: set active bounds (except trust region constraint) as fixed variables
      for (size_t variable_index: this->active_set.bounds.at_lower_bound) {
         const double lb = problem.variable_lower_bound(variable_index);
         if (-INF<double> < lb && std::abs(this->direction_primals[variable_index] + this->trust_region_radius) >= this->activity_tolerance) {
            this->variables_lower_bounds[variable_index] = this->variables_upper_bounds[variable_index] = lb;
         }
      }
      for (size_t variable_index: this->active_set.bounds.at_upper_bound) {
         const double ub = problem.variable_upper_bound(variable_index);
         if (ub < INF<double> && std::abs(this->direction_primals[variable_index] - this->trust_region_radius) >= this->activity_tolerance) {
            this->variables_lower_bounds[variable_index] = this->variables_upper_bounds[variable_index] = ub;
         }
      }

      // constraints: initialize all constraint bounds to {-inf, +inf}
      for (size_t constraint_index: Range(problem.number_constraints)) {
         this->constraints_lower_bounds[constraint_index] = -INF<double>;
         this->constraints_upper_bounds[constraint_index] = INF<double>;
      }
      // constraints: set active bounds as equality constraints
      for (size_t constraint_index: this->active_set.constraints.at_lower_bound) {
         this->constraints_lower_bounds[constraint_index] = this->constraints_upper_bounds[constraint_index] =
            problem.constraint_lower_bound(constraint_index);
      }
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