#include <cmath>
#include "LineSearch.hpp"
#include "Logger.hpp"

LineSearch::LineSearch(GlobalizationStrategy& globalization_strategy, int max_iterations, double ratio):
		GlobalizationMechanism(globalization_strategy, max_iterations), ratio(ratio) {
}

Iterate LineSearch::compute_iterate(Problem& problem, Iterate& current_iterate) {
	this->number_iterations = 0;
	
	/* compute a trial direction */
	LocalSolution solution = this->globalization_strategy.compute_step(problem, current_iterate, INFINITY);

	/* length follows the following sequence: 1, ratio, ratio^2, ratio^3, ... */
	double step_length = 1.;
	bool is_accepted = false;
	while (!this->termination(is_accepted, this->number_iterations)) {
		this->number_iterations++;
		DEBUG << "\n\tLINE SEARCH iteration " << this->number_iterations << ", step_length " << step_length << "\n";
			
		try {
			/* check whether the trial step is accepted */
			is_accepted = this->globalization_strategy.check_step(problem, current_iterate, solution, step_length);
		}
		catch (const std::invalid_argument& e) {
			WARNING << RED << e.what() << RESET "\n";
			is_accepted = false;
		}
		
		if (is_accepted) {
			DEBUG << CYAN "LS trial point accepted\n" RESET;
			/* print summary */
			// TODO handle the case where the solution is a new point: what???
			INFO << "minor: " << std::fixed << this->number_iterations << "\t";
			INFO << "step length: " << step_length << "\t";
			INFO << "step norm: " << std::fixed << (step_length*solution.norm) << "\t";
		}
		else {
			/* decrease the step length */
			step_length *= this->ratio;
		}
	}
	return current_iterate;
}

bool LineSearch::termination(bool is_accepted, int iteration) {
	if (is_accepted) {
		return true;
	}
	else if (this->max_iterations < iteration) {
		throw std::out_of_range("Iteration limit reached");
	}
	/* step length gets too small */
	//if (step_length < 1e-16) {
	//	throw std::out_of_range("Step length became too small");
	//}
	return false;
}
