#include "GlobalizationMechanism.hpp"

GlobalizationMechanism::GlobalizationMechanism(GlobalizationStrategy& globalization_strategy, double tolerance, int max_iterations):
globalization_strategy(globalization_strategy), tolerance(tolerance), max_iterations(max_iterations), number_iterations(0) {
}

GlobalizationMechanism::~GlobalizationMechanism() {
}

std::optional<std::pair<Iterate, Direction> > GlobalizationMechanism::find_first_acceptable_direction_(Statistics& statistics, Problem& problem, Iterate& current_iterate, std::vector<Direction>& directions, double step_length) {
    for (Direction& direction: directions) {
        try {
            std::optional<Iterate> acceptance_check = this->globalization_strategy.check_acceptance(statistics, problem, current_iterate, direction, step_length);
            if (acceptance_check.has_value()) {
                Iterate& accepted_iterate = acceptance_check.value();
                double step_norm = step_length*direction.norm;
                accepted_iterate.status = this->compute_termination_status_(problem, accepted_iterate, step_norm, direction.objective_multiplier);
                /* print summary */
                this->print_acceptance_();
                return std::pair<Iterate, Direction>(accepted_iterate, direction);
            }
        }
        catch (const NumericalError& e) {}
    }
    return std::nullopt;
}

TerminationStatus GlobalizationMechanism::compute_termination_status_(Problem& problem, Iterate& current_iterate, double step_norm, double objective_multiplier) {
    TerminationStatus status = NOT_OPTIMAL;

    if (current_iterate.residuals.complementarity <= this->tolerance * (current_iterate.x.size() + problem.number_constraints)) {
        // feasible and KKT point
        if (current_iterate.residuals.constraints <= this->tolerance * current_iterate.x.size()) {
            if (current_iterate.residuals.KKT <= this->tolerance * std::sqrt(current_iterate.x.size())) {
                status = KKT_POINT;
            }
        }
            // infeasible and FJ point
        else {
            if (current_iterate.residuals.FJ <= this->tolerance * std::sqrt(current_iterate.x.size())) {
                status = FJ_POINT;
            }
        }
    }
    else if (step_norm <= this->tolerance / 100.) {
        if (current_iterate.residuals.constraints <= this->tolerance * current_iterate.x.size()) {
            status = FEASIBLE_SMALL_STEP;
        }
        else {
            status = INFEASIBLE_SMALL_STEP;
        }
    }

    // if convergence, correct the multipliers
    if (status != NOT_OPTIMAL && 0. < objective_multiplier) {
        for (int j = 0; j < problem.number_constraints; j++) {
            current_iterate.multipliers.constraints[j] /= objective_multiplier;
        }
        for (unsigned int i = 0; i < current_iterate.x.size(); i++) {
            current_iterate.multipliers.lower_bounds[i] /= objective_multiplier;
            current_iterate.multipliers.upper_bounds[i] /= objective_multiplier;
        }
    }
    return status;
}