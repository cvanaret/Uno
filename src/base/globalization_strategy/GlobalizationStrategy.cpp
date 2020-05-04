#include <ostream>
#include <cmath>
#include "GlobalizationStrategy.hpp"

GlobalizationStrategy::GlobalizationStrategy(Subproblem& subproblem, double tolerance) : subproblem(subproblem), tolerance(tolerance) {
}

GlobalizationStrategy::~GlobalizationStrategy() {
}

OptimalityStatus GlobalizationStrategy::compute_status(Problem& problem, Iterate& current_iterate, double step_norm, double objective_multiplier) {
    OptimalityStatus status = NOT_OPTIMAL;
    
    if (current_iterate.KKT_residual <= this->tolerance * std::sqrt(current_iterate.x.size()) && current_iterate.complementarity_residual <= this->tolerance * (current_iterate.x.size() + problem.number_constraints)) {
        if (current_iterate.constraint_residual <= this->tolerance * current_iterate.x.size()) {
            status = KKT_POINT;
        }
        else {
             status = FJ_POINT;
        }
    }
    else if (step_norm <= this->tolerance / 100.) {
        if (current_iterate.constraint_residual <= this->tolerance * current_iterate.x.size()) {
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