#include "TrustRegion.hpp"
#include "Utils.hpp"
#include "Logger.hpp"

TrustRegion::TrustRegion(GlobalizationStrategy& globalization_strategy, double initial_radius, int max_iterations):
		GlobalizationMechanism(globalization_strategy, max_iterations), radius(initial_radius) {
	this->radius_max_ = initial_radius;
	this->radius_min_ = initial_radius;
	this->radius_sum_ = initial_radius;
}

Iterate TrustRegion::compute_iterate(Problem& problem, Iterate& current_point) {
	bool is_accepted = false;
	this->number_iterations = 0;
	
	DEBUG << "STARTING TR with radius " << this->radius << "\n";
	/* reset radius (should only increase, if first minor successful) */
	// TODO
	
	while (!this->termination_criterion(is_accepted, this->number_iterations, this->radius)) {
		try {
			this->number_iterations++;
			DEBUG << "\n\tTRUST REGION iteration " << this->number_iterations << ", radius " << this->radius << "\n";
			/* keep track of min/max/average radius */
			this->record_radius(this->radius);
			
			/* compute the step within trust region */
			LocalSolution solution = this->globalization_strategy.compute_step(problem, current_point, this->radius);
			
			/* set multipliers of active trust region to 0 */
			this->correct_multipliers(problem, solution);
			
			/* generate a trial step and check whether it is accepted */
			is_accepted = this->globalization_strategy.check_step(problem, current_point, solution);
			
			if (is_accepted) {
				DEBUG << CYAN "TR trial point accepted\n" RESET;
				/* print summary */
				INFO << "minor: " << std::fixed << this->number_iterations << "\t";
				INFO << "radius: " << std::fixed << this->radius << "\t";
				INFO << "step norm: " << std::fixed << solution.norm << "\t";
				
				/* increase the radius if trust region is active, otherwise keep the same radius */
				if (solution.norm == this->radius) {
					this->radius *= 2.;
				}
			}
			else {
				/* decrease the radius */
				this->radius = std::max(std::min(this->radius, solution.norm)/2., 1e-4);
				
			}
		}
		catch (const std::invalid_argument& e) {
			WARNING << RED << e.what() << "\n" RESET;
			this->radius /= 2.;
		}
	}

	if (this->max_iterations < this->number_iterations) {
		throw std::out_of_range("Iteration limit reached");
	}
	
	return current_point;
}

void TrustRegion::correct_multipliers(Problem& problem, LocalSolution& solution) {
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

bool TrustRegion::termination_criterion(bool is_accepted, int iteration, double radius) {
	/* 1e-16: something like machine precision */
	return is_accepted || this->max_iterations < iteration || radius <= 1e-18;
}

void TrustRegion::record_radius(double radius) {
	this->radius_max_ = std::max(this->radius_max_, radius);
	this->radius_min_ = std::min(this->radius_min_, radius);
	this->radius_sum_ += radius;
	return;
}
