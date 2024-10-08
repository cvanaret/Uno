// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "InequalityConstrainedSubproblem.hpp"
#include "optimization/Iterate.hpp"
#include "reformulation/OptimizationProblem.hpp"

namespace uno {
   void InequalityConstrainedSubproblem::variables_bounds(Vector<double>& variables_lower_bounds, Vector<double>& variables_upper_bounds) const {
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

   void InequalityConstrainedSubproblem::linearized_constraint_bounds(const Vector<double>& constraints, Vector<double>& constraints_lower_bounds,
         Vector<double>& constraints_upper_bounds) const {
      for (size_t constraint_index: Range(this->problem.number_constraints)) {
         constraints_lower_bounds[constraint_index] = problem.constraint_lower_bound(constraint_index) - constraints[constraint_index];
         constraints_upper_bounds[constraint_index] = problem.constraint_upper_bound(constraint_index) - constraints[constraint_index];
      }
   }
} // namespace