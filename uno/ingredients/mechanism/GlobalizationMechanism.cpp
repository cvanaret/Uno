#include "GlobalizationMechanism.hpp"

GlobalizationMechanism::GlobalizationMechanism(ConstraintRelaxationStrategy& constraint_relaxation_strategy, int max_iterations) :
      relaxation_strategy(constraint_relaxation_strategy),
      max_iterations(max_iterations) {
   std::cout << "RELAX: " << this->relaxation_strategy.get_number_variables() << ",  " << this->relaxation_strategy.get_number_constraints() << "\n";
}

Iterate GlobalizationMechanism::assemble_trial_iterate(const Iterate& current_iterate, Direction& direction, double step_length) {
   // TODO do not reevaluate if d = 0
   Iterate trial_iterate(direction.x.size(), direction.multipliers.constraints.size());
   add_vectors(current_iterate.x, direction.x, step_length, trial_iterate.x);
   add_vectors(current_iterate.multipliers.constraints, direction.multipliers.constraints, step_length, trial_iterate.multipliers.constraints);
   copy_from(trial_iterate.multipliers.lower_bounds, direction.multipliers.lower_bounds);
   copy_from(trial_iterate.multipliers.upper_bounds, direction.multipliers.upper_bounds);
   return trial_iterate;
}

int GlobalizationMechanism::get_hessian_evaluation_count() const {
   return this->relaxation_strategy.get_hessian_evaluation_count();
}

int GlobalizationMechanism::get_number_subproblems_solved() const {
   return this->relaxation_strategy.get_number_subproblems_solved();
}

void GlobalizationMechanism::print_warning(const char* message) {
   WARNING << RED << message << RESET << "\n";
}