#include <iostream>
#include <cmath>
#include "TrustLineSearch.hpp"
#include "AMPLModel.hpp"
#include "TrustRegion.hpp"

TrustLineSearch::TrustLineSearch(GlobalizationStrategy& globalization_strategy, double tolerance, double initial_radius, int max_iterations, double ratio) :
GlobalizationMechanism(globalization_strategy, tolerance, max_iterations), ratio(ratio), radius(initial_radius), activity_tolerance_(1e-6) {
}

Iterate TrustLineSearch::initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers) {
    return this->globalization_strategy.initialize(problem, x, multipliers);
}

Iterate TrustLineSearch::compute_acceptable_iterate(Problem& problem, Iterate& current_iterate) {
    bool is_accepted = false;
    this->number_iterations = 0;
    bool linesearch_failed = false;

    while (!this->termination_(is_accepted, this->number_iterations)) {
        try {
            /* compute a trial direction */
            SubproblemSolution solution = this->globalization_strategy.subproblem.compute_step(problem, current_iterate, this->radius);
            
            /* fail if direction is not a descent direction */
            if (0. < dot(solution.x, current_iterate.objective_gradient)) {
                INFO << RED "Trust-line-search direction is not a descent direction\n" RESET;
                linesearch_failed = true;
            }
            else {
                /* set multipliers of active trust region to 0 */
                TrustRegion::correct_active_set(solution, this->radius);

                /* length follows the following sequence: 1, ratio, ratio^2, ratio^3, ... */
                double step_length = 1.;
                while (!this->termination_(is_accepted, this->number_iterations)) {
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
                        DEBUG << CYAN "TLS trial point accepted\n" RESET;
                        current_iterate.status = this->compute_status_(problem, current_iterate, step_length*solution.norm, solution.objective_multiplier);
                        /* print summary */
                        INFO << "minor: " << this->number_iterations << "\t";
                        INFO << "radius: " << this->radius << "\t";
                        INFO << "step length: " << step_length << "\t";
                        INFO << "step norm: " << solution.norm << "\t";

                        /* increase the radius if trust region is active, otherwise keep the same radius */
                        if (solution.norm >= this->radius - this->activity_tolerance_) {
                            this->radius *= 2.;
                        }
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
        }
        catch (const IEEE_Error& e) {
            DEBUG << RED << e.what() << "\n" RESET;
            linesearch_failed = true;
        }
        /* if the line search failed, reduce the trust region radius */
        if (linesearch_failed) {
            if (this->radius == INFINITY) {
                this->radius = 20.;
            }
            this->radius /= 2.;
            this->number_iterations = 0;
        }
    }

    return current_iterate;
}

bool TrustLineSearch::termination_(bool is_accepted, int iteration) {
    if (is_accepted) {
        return true;
    }
    else if (this->max_iterations < iteration) {
        throw std::out_of_range("Trust-line-search iteration limit reached");
    } /* radius gets too small */
    else if (this->radius < 1e-16) { /* 1e-16: something like machine precision */
        throw std::out_of_range("Trust-line-search radius became too small");
    }
    return false;
}
