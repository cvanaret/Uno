// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ActiveSetSubproblem.hpp"

ActiveSetSubproblem::ActiveSetSubproblem(size_t max_number_variables, size_t max_number_constraints):
      Subproblem(max_number_variables, max_number_constraints),
      initial_point(max_number_variables),
      direction_bounds(max_number_variables),
      linearized_constraint_bounds(max_number_constraints) {
}

void ActiveSetSubproblem::generate_initial_iterate(const NonlinearProblem& /*problem*/, Iterate& /*initial_iterate*/) {
}

void ActiveSetSubproblem::set_initial_point(const std::vector<double>& initial_point) {
   copy_from(this->initial_point, initial_point);
}

void ActiveSetSubproblem::initialize_feasibility_problem() {
   // do nothing
}

void ActiveSetSubproblem::set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate) {
   // set the values of the elastic variables
   const auto elastic_setting_function = [&](Iterate& iterate, size_t /*j*/, size_t elastic_index, double /*jacobian_coefficient*/) {
      iterate.primals[elastic_index] = 0.;
      iterate.multipliers.lower_bounds[elastic_index] = 1.;
   };
   problem.set_elastic_variable_values(current_iterate, elastic_setting_function);
}

void ActiveSetSubproblem::exit_feasibility_problem(const NonlinearProblem& /*problem*/, Iterate& /*trial_iterate*/) {
   // do nothing
}

void ActiveSetSubproblem::set_direction_bounds(const NonlinearProblem& problem, const Iterate& current_iterate) {
   // bounds of original variables intersected with trust region
   for (size_t i: Range(problem.get_number_original_variables())) {
      double lb = std::max(-this->trust_region_radius, problem.get_variable_lower_bound(i) - current_iterate.primals[i]);
      double ub = std::min(this->trust_region_radius, problem.get_variable_upper_bound(i) - current_iterate.primals[i]);
      this->direction_bounds[i] = {lb, ub};
   }
   // bounds of additional variables (no trust region!)
   for (size_t i: Range(problem.get_number_original_variables(), problem.number_variables)) {
      const double lb = problem.get_variable_lower_bound(i) - current_iterate.primals[i];
      const double ub = problem.get_variable_upper_bound(i) - current_iterate.primals[i];
      this->direction_bounds[i] = {lb, ub};
   }
}

void ActiveSetSubproblem::set_linearized_constraint_bounds(const NonlinearProblem& problem, const std::vector<double>& current_constraints) {
   for (size_t j: Range(problem.number_constraints)) {
      const double lb = problem.get_constraint_lower_bound(j) - current_constraints[j];
      const double ub = problem.get_constraint_upper_bound(j) - current_constraints[j];
      this->linearized_constraint_bounds[j] = {lb, ub};
   }
}

void ActiveSetSubproblem::compute_dual_displacements(const NonlinearProblem& problem, const Iterate& current_iterate, Direction& direction) {
   // compute dual *displacements* (note: active-set methods usually compute the new duals, not the displacements)
   for (size_t j: Range(problem.number_constraints)) {
      direction.multipliers.constraints[j] -= current_iterate.multipliers.constraints[j];
   }
   for (size_t i: Range(problem.number_variables)) {
      direction.multipliers.lower_bounds[i] -= current_iterate.multipliers.lower_bounds[i];
      direction.multipliers.upper_bounds[i] -= current_iterate.multipliers.upper_bounds[i];
   }
}

void ActiveSetSubproblem::set_auxiliary_measure(const NonlinearProblem& /*problem*/, Iterate& iterate) {
   iterate.progress.auxiliary_terms = 0.;
}

double ActiveSetSubproblem::generate_predicted_auxiliary_reduction_model(const NonlinearProblem& /*problem*/,
      const Iterate& /*current_iterate*/, const Direction& /*direction*/, double /*step_length*/) const {
   return 0.;
   //}, "0"};
}

void ActiveSetSubproblem::postprocess_iterate(const NonlinearProblem& /*problem*/, Iterate& /*iterate*/) {
}