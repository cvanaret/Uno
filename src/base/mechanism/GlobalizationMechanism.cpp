#include "GlobalizationMechanism.hpp"

GlobalizationMechanism::GlobalizationMechanism(ConstraintRelaxationStrategy& constraint_relaxation_strategy, int max_iterations) :
      relaxation_strategy(constraint_relaxation_strategy), max_iterations(max_iterations), number_iterations(0) {
}

Iterate GlobalizationMechanism::assemble_trial_iterate(const Iterate& current_iterate, Direction& direction, double step_length) {
   // TODO do not reevaluate if ||d|| = 0
   add_vectors(current_iterate.x, direction.x, step_length, this->trial_primals_);
   add_vectors(current_iterate.multipliers.constraints, direction.multipliers.constraints, step_length, this->trial_duals_);
   direction.multipliers.constraints = this->trial_duals_;
   Iterate trial_iterate(this->trial_primals_, direction.multipliers);
   return trial_iterate;
}

int GlobalizationMechanism::get_hessian_evaluation_count() const {
   return this->relaxation_strategy.get_hessian_evaluation_count();
}

int GlobalizationMechanism::get_number_subproblems_solved() const {
   return this->relaxation_strategy.get_number_subproblems_solved();
}

void GlobalizationMechanism::print_acceptance_(const Iterate& iterate) {
   DEBUG << CYAN "trial point accepted\n" RESET;
   DEBUG << "Residuals: ||c|| = " << iterate.errors.constraints << ", KKT = " << iterate.errors.KKT
         << ", complementarity = " << iterate.errors.complementarity << "\n";
}

void GlobalizationMechanism::print_warning_(const char* message) {
   WARNING << RED << message << RESET << "\n";
}