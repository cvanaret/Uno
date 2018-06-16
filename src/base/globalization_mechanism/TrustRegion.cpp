#include <cmath>
#include "TrustRegion.hpp"
#include "Utils.hpp"
#include "Logger.hpp"

TrustRegion::TrustRegion(GlobalizationStrategy& globalization_strategy, double initial_radius, int max_iterations):
		GlobalizationMechanism(globalization_strategy, max_iterations), radius(initial_radius) {
	this->radius_max_ = -INFINITY;
	this->radius_min_ = INFINITY;
	this->radius_sum_ = 0;
	
	this->activity_tolerance_ = 1e-6;
}

Iterate TrustRegion::compute_iterate(Problem& problem, Iterate& current_iterate) {
	bool is_accepted = false;
	this->number_iterations = 0;
	
	while (!this->termination(is_accepted, this->number_iterations, this->radius)) {
		try {
			this->number_iterations++;
			/* keep track of min/max/average radius */
			this->record_radius(this->radius);
			
			DEBUG << "\n\tTRUST REGION iteration " << this->number_iterations << ", radius " << this->radius << "\n";
			
			/* compute the step within trust region */
			LocalSolution solution = this->globalization_strategy.compute_step(problem, current_iterate, this->radius);
			
			/* set multipliers of active trust region to 0 */
			this->correct_multipliers(problem, solution);
			
			/* check whether the trial step is accepted */
			is_accepted = this->globalization_strategy.check_step(problem, current_iterate, solution);
			
			if (is_accepted) {
				DEBUG << CYAN "TR trial point accepted\n" RESET;
				/* print summary */
				INFO << "minor: " << std::fixed << this->number_iterations << "\t";
				INFO << "radius: " << std::fixed << this->radius << "\t";
				INFO << "step norm: " << std::fixed << solution.norm << "\t";
				
				/* increase the radius if trust region is active, otherwise keep the same radius */
				if (solution.norm >= this->radius - this->activity_tolerance_) {
					this->radius *= 2.;
				}
			}
			else {
				/* decrease the radius */
				this->radius = std::min(this->radius, solution.norm)/2.;
			}
		}
		catch (const std::invalid_argument& e) {
			WARNING << RED << e.what() << RESET "\n";
			this->radius /= 2.;
		}
	}
	return current_iterate;
}

void TrustRegion::correct_multipliers(Problem& problem, LocalSolution& solution) {
	/* set multipliers for bound constraints active at trust region to 0 */
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

bool TrustRegion::termination(bool is_accepted, int iteration, double radius) {
	if (is_accepted) {
		return true;
	}
	else if (this->max_iterations < iteration) {
		throw std::out_of_range("Trust-region iteration limit reached");
	}
	/* radius gets too small */
	else if (this->radius < 1e-16) { /* 1e-16: something like machine precision */
		throw std::out_of_range("Trust-region radius became too small");
	}
	return false;
}

void TrustRegion::record_radius(double radius) {
	this->radius_max_ = std::max(this->radius_max_, radius);
	this->radius_min_ = std::min(this->radius_min_, radius);
	this->radius_sum_ += radius;
	return;
}
