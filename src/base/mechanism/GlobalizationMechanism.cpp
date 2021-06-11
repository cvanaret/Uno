#include "GlobalizationMechanism.hpp"

GlobalizationMechanism::GlobalizationMechanism(ConstraintRelaxationStrategy& constraint_relaxation_strategy, int max_iterations) :
      relaxation_strategy(constraint_relaxation_strategy), max_iterations(max_iterations), number_iterations(0) {
}

void GlobalizationMechanism::print_acceptance_(const Iterate& iterate) {
   DEBUG << CYAN "trial point accepted\n" RESET;
   DEBUG << "Residuals: ||c|| = " << iterate.residuals.constraints << ", KKT = " << iterate.residuals.KKT
         << ", complementarity = " << iterate.residuals.complementarity << "\n";
}

void GlobalizationMechanism::print_warning_(const char* message) {
   WARNING << RED << message << RESET << "\n";
}