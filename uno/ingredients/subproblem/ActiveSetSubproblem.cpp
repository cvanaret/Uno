// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ActiveSetSubproblem.hpp"

ActiveSetSubproblem::ActiveSetSubproblem(size_t max_number_variables, size_t max_number_constraints):
      Subproblem(max_number_variables, max_number_constraints),
      initial_point(max_number_variables),
      variable_displacement_bounds(max_number_variables),
      linearized_constraint_bounds(max_number_constraints) {
}

void ActiveSetSubproblem::initialize(Statistics& /*statistics*/, const NonlinearProblem& /*problem*/, Iterate& /*first_iterate*/) {
}

void ActiveSetSubproblem::set_initial_point(const std::vector<double>& initial_point) {
   copy_from(this->initial_point, initial_point);
}

void ActiveSetSubproblem::prepare_for_feasibility_problem(const NonlinearProblem& /*problem*/, Iterate& /*current_iterate*/) {
   // do nothing
}

void ActiveSetSubproblem::set_elastic_variables(const l1RelaxedProblem& problem, Iterate& current_iterate) {
   // reset the values of the elastic variables
   const auto elastic_setting_function = [&](Iterate& iterate, size_t /*j*/, size_t elastic_index, double /*jacobian_coefficient*/,
         double /*constraint_violation_coefficient*/) {
      iterate.primals[elastic_index] = 0.;
      iterate.multipliers.lower_bounds[elastic_index] = 1.;
   };
   problem.set_elastic_variables(current_iterate, elastic_setting_function);
}

void ActiveSetSubproblem::set_variable_displacement_bounds(const NonlinearProblem& problem, const Iterate& current_iterate) {
   for (size_t i = 0; i < problem.number_variables; i++) {
      const double lb = this->variable_bounds[i].lb - current_iterate.primals[i];
      const double ub = this->variable_bounds[i].ub - current_iterate.primals[i];
      this->variable_displacement_bounds[i] = {lb, ub};
   }
}

void ActiveSetSubproblem::set_linearized_constraint_bounds(const NonlinearProblem& problem, const std::vector<double>& current_constraints) {
   for (size_t j = 0; j < problem.number_constraints; j++) {
      const double lb = problem.get_constraint_lower_bound(j) - current_constraints[j];
      const double ub = problem.get_constraint_upper_bound(j) - current_constraints[j];
      this->linearized_constraint_bounds[j] = {lb, ub};
   }
}

void ActiveSetSubproblem::shift_linearized_constraint_bounds(const NonlinearProblem& problem, const std::vector<double>& trial_constraints) {
   // shift the RHS with the values of the constraints at the trial iterate
   for (size_t j = 0; j < problem.number_constraints; j++) {
      this->linearized_constraint_bounds[j].lb -= trial_constraints[j];
      this->linearized_constraint_bounds[j].ub -= trial_constraints[j];
   }
}

void ActiveSetSubproblem::compute_dual_displacements(const NonlinearProblem& problem, const Iterate& current_iterate, Direction& direction) {
   // compute dual *displacements* (note: active-set methods usually compute the new duals, not the displacements)
   for (size_t j = 0; j < problem.number_constraints; j++) {
      direction.multipliers.constraints[j] -= current_iterate.multipliers.constraints[j];
   }
}

void ActiveSetSubproblem::set_optimality_measure(const NonlinearProblem& problem, Iterate& iterate) {
   // optimality measure: original objective value
   iterate.evaluate_objective(problem.model);
   iterate.nonlinear_progress.optimality = iterate.original_evaluations.objective;
}

void ActiveSetSubproblem::postprocess_accepted_iterate(const NonlinearProblem& /*problem*/, Iterate& /*iterate*/) {
}