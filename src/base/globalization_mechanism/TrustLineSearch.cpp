#include <iostream>
#include <cmath>
#include "TrustLineSearch.hpp"

TrustLineSearch::TrustLineSearch(GlobalizationStrategy& globalization_strategy, double initial_radius, int max_iterations, double ratio):
		GlobalizationMechanism(globalization_strategy, max_iterations), ratio(ratio), radius(initial_radius) {
			
	this->activity_tolerance_ = 1e-6;
}

Iterate TrustLineSearch::compute_iterate(Problem& problem, Iterate& current_iterate) {
	this->number_iterations = 0;

	bool is_accepted = false;
	bool linesearch_failed = false;
	while (!this->termination_criterion(is_accepted, this->number_iterations)) {
		try {
			if (this->radius == 0.) {
				throw std::out_of_range("radius is zero");
			}
			/* compute a trial direction */
			LocalSolution solution = this->globalization_strategy.compute_step(problem, current_iterate, this->radius);
			
			/* set multipliers of active trust region to 0 */
			this->correct_multipliers(problem, solution);

			/* length follows the following sequence: 1, ratio, ratio^2, ratio^3, ... */
			double step_length = 1.;
			while (!this->termination_criterion(is_accepted, this->number_iterations)) {
				this->number_iterations++;
				DEBUG << "\n\tTRUST LINE SEARCH iteration " << this->number_iterations << ", radius " << this->radius << ", step_length " << step_length << "\n";
				
				try {
					/* check whether the trial step is accepted */
					is_accepted = this->globalization_strategy.check_step(problem, current_iterate, solution, step_length);
				}
				catch (const std::invalid_argument& e) {
					DEBUG << RED << e.what() << "\n" RESET;
					is_accepted = false;
				}
			
				if (is_accepted) {
					DEBUG << "LS trial point accepted\n";
					
					/* increase the radius if trust region is active, otherwise keep the same radius */
					// TODO handle the case where the solution is the new point
					//if (solution.norm == this->radius - this->activity_tolerance_) {
					//	std::cout << this->radius << " RADIUS MULTIPLIED BY 2\n";
					//	this->radius *= 2.;
					//}
				}
				else {
					/* decrease alpha */
					step_length *= this->ratio;
				}
			}
			if (this->max_iterations < this->number_iterations) {
				linesearch_failed = true;
			}
		}
		catch (const std::invalid_argument& e) {
			DEBUG << RED << e.what() << "\n" RESET;
			linesearch_failed = true;
		}
		/* if the line search failed, reduce the trust region radius */
		if (linesearch_failed) {
			std::cout << this->radius << " RADIUS DIVIDED BY 2\n";
			this->radius /= 2.;
			this->number_iterations = 0;
		}
		
		if (this->radius < 1e-4) {
			this->radius = 1e-4;
		}
	}

	return current_iterate;
}

void TrustLineSearch::correct_multipliers(Problem& problem, LocalSolution& solution) {
	for (unsigned int k = 0; k < solution.active_set.at_upper_bound.size(); k++) {
		int i = solution.active_set.at_upper_bound[k];
		if (i < problem.number_variables && solution.x[i] == this->radius) {
			solution.multipliers[i] = 0.;
		}
	}
	for (unsigned int k = 0; k < solution.active_set.at_lower_bound.size(); k++) {
		int i = solution.active_set.at_lower_bound[k];
		if (i < problem.number_variables && solution.x[i] == -this->radius) {
			solution.multipliers[i] = 0.;
		}
	}
	return;
}

bool TrustLineSearch::termination_criterion(bool is_accepted, int iteration) {
	return is_accepted || this->max_iterations < iteration;
}
