#include <cassert>
#include "Subproblem.hpp"
#include "linear_algebra/SparseVector.hpp"

Subproblem::Subproblem(const NonlinearProblem& problem, SecondOrderCorrection soc_strategy):
      objective_gradient(problem.number_variables), constraints(problem.number_constraints), constraint_jacobian(problem.number_constraints),
      soc_strategy(soc_strategy), variable_bounds(problem.number_variables),
      direction(problem.number_variables, problem.number_constraints) {
   for (auto& constraint_gradient: this->constraint_jacobian) {
      constraint_gradient.reserve(problem.number_variables);
   }
   // register the variables bounds
   for (size_t i = 0; i < problem.number_variables; i++) {
      this->variable_bounds[i] = {problem.get_variable_lower_bound(i), problem.get_variable_upper_bound(i)};
   }
}

void Subproblem::set_variable_bounds(const NonlinearProblem& problem, const Iterate& current_iterate, double trust_region_radius) {
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