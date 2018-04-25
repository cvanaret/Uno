#include <cmath>
#include "LineSearch.hpp"
#include "Logger.hpp"

LineSearch::LineSearch(GlobalizationStrategy& globalization_strategy, int max_iterations):
		GlobalizationMechanism(globalization_strategy, max_iterations) {
}

// TODO catch exceptions
Iterate LineSearch::compute_iterate(Problem& problem, Iterate& current_point) {
	this->number_iterations = 0;
	double ratio = 0.5; // in ]0, 1[

	/* compute a trial direction */
	LocalSolution solution = this->globalization_strategy.compute_step(problem, current_point, INFINITY);

	/* length follows the following sequence: 1, ratio, ratio^2, ratio^3, ... */
	double step_length = 1.;
	bool is_accepted = false;
	while (!this->termination_criterion(is_accepted, this->number_iterations)) {
		this->number_iterations++;
		DEBUG << "\n\tLINE SEARCH iteration " << this->number_iterations << ", step_length " << step_length << "\n";
			
		try {
			/* check whether the trial step is accepted */
			is_accepted = this->globalization_strategy.check_step(problem, current_point, solution, step_length);
		}
		catch (const std::invalid_argument& e) {
			//WARNING << RED << e.what() << "\n" RESET;
			is_accepted = false;
		}
		
		if (is_accepted) {
			DEBUG << CYAN "LS trial point accepted\n" RESET;
		}
		else {
			/* decrease alpha */
			step_length *= ratio;
			if (step_length < 1e-8) {
				step_length = 0.;
			}
		}
	}

	if (this->max_iterations < this->number_iterations) {
		throw std::out_of_range("Iteration limit reached");
	}
	
	/* print summary */
	// TODO handle the case where the solution is a new point
	double step_norm = solution.norm*step_length;
	INFO << "minor: " << std::fixed << this->number_iterations << "\t";
	INFO << "step length: " << std::fixed << step_length << "\t";
	INFO << "step norm: " << std::fixed << step_norm << "\t";

	return current_point;
}

bool LineSearch::termination_criterion(bool is_accepted, int iteration) {
	return is_accepted || this->max_iterations < iteration;
}
