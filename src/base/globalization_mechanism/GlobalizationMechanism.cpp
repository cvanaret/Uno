#include "GlobalizationMechanism.hpp"

GlobalizationMechanism::GlobalizationMechanism(GlobalizationStrategy& globalization_strategy, int max_iterations):
		globalization_strategy(globalization_strategy), max_iterations(max_iterations), number_iterations(0) {
}

GlobalizationMechanism::~GlobalizationMechanism() {
}

void GlobalizationMechanism::initialize(Problem& problem, Iterate& current_iterate) {
	this->globalization_strategy.initialize(problem, current_iterate);
	return;
}
