// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "InequalityConstrainedMethod.hpp"

InequalityConstrainedMethod::InequalityConstrainedMethod(size_t max_number_variables, size_t max_number_constraints):
      Subproblem(max_number_variables, max_number_constraints),
      initial_point(max_number_variables),
      direction_bounds(max_number_variables),
      linearized_constraint_bounds(max_number_constraints) {
}

void InequalityConstrainedMethod::generate_initial_iterate(const NonlinearProblem& /*problem*/, Iterate& /*initial_iterate*/) {
}

void InequalityConstrainedMethod::set_initial_point(const std::vector<double>& initial_point) {
   copy_from(this->initial_point, initial_point);
}

void InequalityConstrainedMethod::initialize_feasibility_problem() {
   // do nothing
}

void InequalityConstrainedMethod::set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate) {
   // set the values of the elastic variables
   const auto elastic_setting_function = [&](Iterate& iterate, size_t /*j*/, size_t elastic_index, double /*jacobian_coefficient*/) {
      iterate.primals[elastic_index] = 0.;
      iterate.multipliers.lower_bounds[elastic_index] = 1.;
   };
   problem.set_elastic_variable_values(current_iterate, elastic_setting_function);
}

void InequalityConstrainedMethod::exit_feasibility_problem(const NonlinearProblem& /*problem*/, Iterate& /*trial_iterate*/) {
   // do nothing
}

void InequalityConstrainedMethod::set_direction_bounds(const NonlinearProblem& problem, const Iterate& current_iterate) {
   // bounds of original variables intersected with trust region
   for (size_t variable_index: Range(problem.get_number_original_variables())) {
      double lb = std::max(-this->trust_region_radius, problem.get_variable_lower_bound(variable_index) - current_iterate.primals[variable_index]);
      double ub = std::min(this->trust_region_radius, problem.get_variable_upper_bound(variable_index) - current_iterate.primals[variable_index]);
      this->direction_bounds[variable_index] = {lb, ub};
   }
   // bounds of additional variables (no trust region!)
   for (size_t variable_index: Range(problem.get_number_original_variables(), problem.number_variables)) {
      const double lb = problem.get_variable_lower_bound(variable_index) - current_iterate.primals[variable_index];
      const double ub = problem.get_variable_upper_bound(variable_index) - current_iterate.primals[variable_index];
      this->direction_bounds[variable_index] = {lb, ub};
   }
}

void InequalityConstrainedMethod::set_linearized_constraint_bounds(const NonlinearProblem& problem, const std::vector<double>& current_constraints) {
   for (size_t constraint_index: Range(problem.number_constraints)) {
      const double lb = problem.get_constraint_lower_bound(constraint_index) - current_constraints[constraint_index];
      const double ub = problem.get_constraint_upper_bound(constraint_index) - current_constraints[constraint_index];
      this->linearized_constraint_bounds[constraint_index] = {lb, ub};
   }
}

void InequalityConstrainedMethod::compute_dual_displacements(const NonlinearProblem& problem, const Iterate& current_iterate, Direction& direction) {
   // compute dual *displacements* (note: active-set methods usually compute the new duals, not the displacements)
   for (size_t constraint_index: Range(problem.number_constraints)) {
      direction.multipliers.constraints[constraint_index] -= current_iterate.multipliers.constraints[constraint_index];
   }
   for (size_t variable_index: Range(problem.number_variables)) {
      direction.multipliers.lower_bounds[variable_index] -= current_iterate.multipliers.lower_bounds[variable_index];
      direction.multipliers.upper_bounds[variable_index] -= current_iterate.multipliers.upper_bounds[variable_index];
   }
}

void InequalityConstrainedMethod::set_auxiliary_measure(const NonlinearProblem& /*problem*/, Iterate& iterate) {
   iterate.progress.auxiliary_terms = 0.;
}

double InequalityConstrainedMethod::compute_predicted_auxiliary_reduction_model(const NonlinearProblem& /*problem*/,
      const Iterate& /*current_iterate*/, const Direction& /*direction*/, double /*step_length*/) const {
   return 0.;
   //}, "0"};
}

void InequalityConstrainedMethod::postprocess_iterate(const NonlinearProblem& /*problem*/, Iterate& /*iterate*/) {
}