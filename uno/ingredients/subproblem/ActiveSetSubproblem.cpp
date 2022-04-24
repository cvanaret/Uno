#include "ActiveSetSubproblem.hpp"

ActiveSetSubproblem::ActiveSetSubproblem(const ReformulatedProblem& problem, SecondOrderCorrection soc_strategy):
      Subproblem(problem, soc_strategy),
      initial_point(problem.number_variables),
      variable_displacement_bounds(problem.number_variables),
      linearized_constraint_bounds(problem.number_constraints) {
}

void ActiveSetSubproblem::initialize(Statistics& /*statistics*/, const ReformulatedProblem& /*problem*/, Iterate& /*first_iterate*/) {
}

void ActiveSetSubproblem::set_initial_point(const std::optional<std::vector<double>>& optional_initial_point) {
   // if provided, set the initial point
   if (optional_initial_point.has_value()) {
      const std::vector<double>& point = optional_initial_point.value();
      copy_from(this->initial_point, point);
   }
   else {
      // otherwise, reset the initial point
      initialize_vector(this->initial_point, 0.);
   }
}

void ActiveSetSubproblem::set_elastic_variables(const l1RelaxedProblem& problem, Iterate& current_iterate) {
   // reset (set to 0) the values of the elastic variables
   problem.set_elastic_variables(current_iterate);
}

void ActiveSetSubproblem::set_variable_displacement_bounds(const ReformulatedProblem& problem, const Iterate& current_iterate) {
   for (size_t i = 0; i < problem.number_variables; i++) {
      const double lb = this->variable_bounds[i].lb - current_iterate.primals[i];
      const double ub = this->variable_bounds[i].ub - current_iterate.primals[i];
      this->variable_displacement_bounds[i] = {lb, ub};
   }
}

void ActiveSetSubproblem::set_linearized_constraint_bounds(const ReformulatedProblem& problem, const std::vector<double>& current_constraints) {
   for (size_t j = 0; j < problem.number_constraints; j++) {
      const double lb = problem.get_constraint_lower_bound(j) - current_constraints[j];
      const double ub = problem.get_constraint_upper_bound(j) - current_constraints[j];
      this->linearized_constraint_bounds[j] = {lb, ub};
   }
}

void ActiveSetSubproblem::compute_dual_displacements(const ReformulatedProblem& problem, const Iterate& current_iterate, Direction& direction) {
   // compute dual *displacements* (note: active-set methods usually compute the new duals, not the displacements)
   for (size_t j = 0; j < problem.number_constraints; j++) {
      direction.multipliers.constraints[j] -= current_iterate.multipliers.constraints[j];
   }
}

double ActiveSetSubproblem::compute_optimality_measure(const ReformulatedProblem& problem, Iterate& iterate) {
   // optimality measure: objective value
   return problem.evaluate_objective(iterate);
}

Direction ActiveSetSubproblem::compute_second_order_correction(const ReformulatedProblem& /*problem*/, Iterate& /*trial_iterate*/) {
   assert(false && "ActiveSetSubproblem::compute_second_order_correction");
}

void ActiveSetSubproblem::postprocess_accepted_iterate(const ReformulatedProblem& /*problem*/, Iterate& /*iterate*/) {
}