#include <cmath>
#include <assert.h>
#include "TrustRegion.hpp"
#include "Vector.hpp"
#include "Logger.hpp"

TrustRegion::TrustRegion(GlobalizationStrategy& globalization_strategy, double tolerance, double initial_radius, int max_iterations):
GlobalizationMechanism(globalization_strategy, tolerance, max_iterations), radius(initial_radius), activity_tolerance_(1e-6) {
}

Iterate TrustRegion::initialize(Statistics& statistics, Problem& problem, std::vector<double>& x, Multipliers& multipliers) {
    statistics.add_column("minor", Statistics::int_width, 2);
    statistics.add_column("TR radius", Statistics::double_width, 30);
    statistics.add_column("step norm", Statistics::double_width, 31);
    return this->globalization_strategy.initialize(statistics, problem, x, multipliers);
}

Iterate TrustRegion::compute_acceptable_iterate(Statistics& statistics, Problem& problem, Iterate& current_iterate) {
    bool is_accepted = false;
    this->number_iterations = 0;

    while (!this->termination_(is_accepted)) {
        try {
            assert (0 < this->radius);
            this->number_iterations++;
            this->print_iteration_();

            /* compute the step within trust region */
            std::vector<Direction> directions = this->globalization_strategy.subproblem.compute_directions(problem, current_iterate, this->radius);
            /* set bound multipliers of active trust region to 0 */
            for (Direction& direction: directions) {
                this->correct_active_set(direction, this->radius);
            }
            
            /* check whether the trial step is accepted */
            std::optional<std::pair<Iterate, double> > acceptance_check = this->find_first_acceptable_direction_(statistics, problem, current_iterate, directions, 1.);
            if (acceptance_check.has_value()) {
                is_accepted = true;
                current_iterate = acceptance_check.value().first;
                double step_norm = acceptance_check.value().second;
                statistics.add_statistic("minor", this->number_iterations);
                statistics.add_statistic("TR radius", this->radius);
                statistics.add_statistic("step norm", step_norm);
                /* increase the radius if trust region is active, otherwise keep the same radius */
                if (step_norm >= this->radius - this->activity_tolerance_) {
                    this->radius *= 2.;
                }
            }
            else {
                /* if the step is rejected, decrease the radius */
                double min_norm = INFINITY;
                for (const Direction& direction: directions) {
                    min_norm = std::min(min_norm, direction.norm);
                }
                this->radius = std::min(this->radius, min_norm) / 2.;
            }
        }
        catch (const IEEE_Error& e) {
            this->print_warning_(e.what());
            /* if an evaluation error occurs, decrease the radius */
            this->radius /= 2.;
        }
    }
    return current_iterate;
}

void TrustRegion::correct_active_set(Direction& direction, double radius) {
    assert (0 < radius);
    /* update active set and set multipliers for bound constraints active at trust region to 0 */
    for (std::set<int>::iterator it = direction.active_set.bounds.at_lower_bound.begin(); it != direction.active_set.bounds.at_lower_bound.end();) {
        int i = *it;
        if (direction.x[i] == -radius) {
            it = direction.active_set.bounds.at_lower_bound.erase(it);
            direction.multipliers.lower_bounds[i] = 0.;
        }
        else {
            ++it;
        }
    }
    for (std::set<int>::iterator it = direction.active_set.bounds.at_upper_bound.begin(); it != direction.active_set.bounds.at_upper_bound.end();) {
        int i = *it;
        if (direction.x[i] == radius) {
            it = direction.active_set.bounds.at_upper_bound.erase(it);
            direction.multipliers.upper_bounds[i] = 0.;
        }
        else {
            ++it;
        }
    }
    return;
}

bool TrustRegion::termination_(bool is_accepted) {
    if (is_accepted) {
        return true;
    }
    else if (this->max_iterations < this->number_iterations) {
        throw std::runtime_error("Trust-region iteration limit reached");
    } /* radius gets too small */
    else if (this->radius < 1e-16) { /* 1e-16: something like machine precision */
        throw std::runtime_error("Trust-region radius became too small");
    }
    return false;
}

void TrustRegion::print_iteration_() {
    DEBUG << "\n\tTRUST REGION iteration " << this->number_iterations << ", radius " << this->radius << "\n";
    return;
}

void TrustRegion::print_acceptance_() {
    DEBUG << CYAN "TR trial point accepted\n" RESET;
    return;
}

void TrustRegion::print_warning_(const char* message) {
    WARNING << RED << message << RESET << "\n";
    return;
}
