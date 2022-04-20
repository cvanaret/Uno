#include <cassert>
#include "Subproblem.hpp"
#include "linear_algebra/SparseVector.hpp"

Subproblem::Subproblem(size_t max_number_variables, size_t number_constraints, SecondOrderCorrection soc_strategy):
      soc_strategy(soc_strategy), variable_bounds(max_number_variables),
      direction(max_number_variables, number_constraints) {
}

void Subproblem::set_variable_bounds(const NonlinearProblem& problem, const Iterate& current_iterate, double trust_region_radius) {
   // bounds intersected with trust region
   // very important: apply the trust region only on the original variables
   for (size_t i = 0; i < problem.get_number_original_variables(); i++) {
      double lb = std::max(current_iterate.x[i] - trust_region_radius, problem.get_variable_lower_bound(i));
      double ub = std::min(current_iterate.x[i] + trust_region_radius, problem.get_variable_upper_bound(i));
      this->variable_bounds[i] = {lb, ub};
   }
   for (size_t i = problem.get_number_original_variables(); i < problem.number_variables; i++) {
      const double lb = problem.get_variable_lower_bound(i);
      const double ub = problem.get_variable_upper_bound(i);
      this->variable_bounds[i] = {lb, ub};
   }
}