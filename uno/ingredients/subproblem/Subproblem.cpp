// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project root for details.

#include <cassert>
#include "Subproblem.hpp"
#include "linear_algebra/SparseVector.hpp"

Subproblem::Subproblem(size_t max_number_variables, size_t max_number_constraints, SecondOrderCorrection soc_strategy):
      objective_gradient(max_number_variables), constraints(max_number_constraints), constraint_jacobian(max_number_constraints),
      soc_strategy(soc_strategy), variable_bounds(max_number_variables),
      direction(max_number_variables, max_number_constraints) {
   for (auto& constraint_gradient: this->constraint_jacobian) {
      constraint_gradient.reserve(max_number_variables);
   }
}

void Subproblem::set_variable_bounds(const ReformulatedProblem& problem, const Iterate& current_iterate, double trust_region_radius) {
   // bounds intersected with trust region
   // very important: apply the trust region only on the original variables
   for (size_t i = 0; i < problem.get_number_original_variables(); i++) {
      double lb = std::max(current_iterate.primals[i] - trust_region_radius, problem.get_variable_lower_bound(i));
      double ub = std::min(current_iterate.primals[i] + trust_region_radius, problem.get_variable_upper_bound(i));
      this->variable_bounds[i] = {lb, ub};
   }
   for (size_t i = problem.get_number_original_variables(); i < problem.number_variables; i++) {
      const double lb = problem.get_variable_lower_bound(i);
      const double ub = problem.get_variable_upper_bound(i);
      this->variable_bounds[i] = {lb, ub};
   }
}