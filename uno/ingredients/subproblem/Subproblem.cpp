// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "Subproblem.hpp"

Subproblem::Subproblem(size_t max_number_variables, size_t max_number_constraints):
      direction(max_number_variables, max_number_constraints),
      evaluations(max_number_variables, max_number_constraints),
      variable_bounds(max_number_variables) {
}

void Subproblem::set_trust_region_radius(double new_trust_region_radius) {
   assert(0. < new_trust_region_radius && "The trust-region radius should be positive.");
   this->trust_region_radius = new_trust_region_radius;
}

void Subproblem::set_variable_bounds(const NonlinearProblem& problem, const Iterate& current_iterate) {
   // bounds intersected with trust region
   // very important: apply the trust region only on the original variables
   for (size_t i: Range(problem.get_number_original_variables())) {
      double lb = std::max(current_iterate.primals[i] - this->trust_region_radius, problem.get_variable_lower_bound(i));
      double ub = std::min(current_iterate.primals[i] + this->trust_region_radius, problem.get_variable_upper_bound(i));
      this->variable_bounds[i] = {lb, ub};
   }
   for (size_t i: Range(problem.get_number_original_variables(), problem.number_variables)) {
      const double lb = problem.get_variable_lower_bound(i);
      const double ub = problem.get_variable_upper_bound(i);
      this->variable_bounds[i] = {lb, ub};
   }
}

void Subproblem::check_unboundedness(const Direction& direction) {
   assert(direction.status != SubproblemStatus::UNBOUNDED_PROBLEM && "The subproblem is unbounded, this should not happen");
}