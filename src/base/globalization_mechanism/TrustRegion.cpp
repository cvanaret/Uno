#include <cmath>
#include "TrustRegion.hpp"
#include "Utils.hpp"
#include "Logger.hpp"

TrustRegion::TrustRegion(GlobalizationStrategy& globalization_strategy, double initial_radius, int max_iterations):
		GlobalizationMechanism(globalization_strategy, max_iterations), radius(initial_radius), activity_tolerance_(1e-6) {
}

void TrustRegion::initialize(Problem& problem, Iterate& current_iterate) {
    this->globalization_strategy.initialize(problem, current_iterate, true);
    return;
}

Iterate TrustRegion::compute_iterate(Problem& problem, Iterate& current_iterate) {
	bool is_accepted = false;
	this->number_iterations = 0;

    while (!this->termination(is_accepted)) {
		try {
            this->number_iterations++;
            this->print_iteration();
			
			/* compute the step within trust region */
			LocalSolution solution = this->globalization_strategy.compute_step(problem, current_iterate, this->radius);

            /* set bound multipliers of active trust region to 0 */
			this->correct_multipliers(problem, solution);
			
			/* check whether the trial step is accepted */
			is_accepted = this->globalization_strategy.check_step(problem, current_iterate, solution);

            if (is_accepted) {
                /* print summary */
                this->print_acceptance(solution.norm);
                
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
            this->print_warning(e.what());
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
            solution.bound_multipliers[i] = 0.;
		}
	}
	for (unsigned int k = 0; k < solution.active_set.at_lower_bound.size(); k++) {
		int i = solution.active_set.at_lower_bound[k];
        if (i < problem.number_variables && solution.x[i] == -this->radius) {
            solution.bound_multipliers[i] = 0.;
		}
	}
	return;
}

bool TrustRegion::termination(bool is_accepted) {
	if (is_accepted) {
		return true;
    }
    else if (this->max_iterations < this->number_iterations) {
		throw std::out_of_range("Trust-region iteration limit reached");
	}
	/* radius gets too small */
	else if (this->radius < 1e-16) { /* 1e-16: something like machine precision */
		throw std::out_of_range("Trust-region radius became too small");
	}
	return false;
}

void TrustRegion::print_iteration() {
    DEBUG << "\n\tTRUST REGION iteration " << this->number_iterations << ", radius " << this->radius << "\n";
    return;
}

void TrustRegion::print_acceptance(double solution_norm) {
    DEBUG << CYAN "TR trial point accepted\n" RESET;
    INFO << "minor: " << this->number_iterations << "\t";
    INFO << "radius: " << this->radius << "\t";
    INFO << "step norm: " << solution_norm << "\t";
    return;
}

void TrustRegion::print_warning(const char* message) {
    WARNING << RED << message << RESET "\n";
    return;
}
