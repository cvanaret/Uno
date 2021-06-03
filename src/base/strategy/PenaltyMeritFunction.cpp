#include <cmath>
#include "Uno.hpp"
#include "PenaltyMeritFunction.hpp"
#include "Logger.hpp"

/*
 * Infeasibility detection and SQP methods for nonlinear optimization 
 * http://epubs.siam.org/doi/pdf/10.1137/080738222
 */

PenaltyMeritFunction::PenaltyMeritFunction(Subproblem& subproblem): GlobalizationStrategy(subproblem), decrease_fraction_(1e-8) {
}

Iterate PenaltyMeritFunction::initialize(Statistics& statistics, Problem& problem, std::vector<double>& x, Multipliers& multipliers) {
    statistics.add_column("penalty param.", Statistics::double_width, 4);
    /* initialize the subproblem */
    Iterate first_iterate = this->subproblem.evaluate_initial_point(problem, x, multipliers);
    this->subproblem.compute_residuals(problem, first_iterate, first_iterate.multipliers, 1.);
    return first_iterate;
}

std::optional<Iterate> PenaltyMeritFunction::check_acceptance(Statistics& statistics, Problem& problem, Iterate& current_iterate, Direction& direction, double step_length) {
    /* check if subproblem definition changed */
    if (this->subproblem.subproblem_definition_changed) {
        this->subproblem.subproblem_definition_changed = false;
        this->subproblem.compute_optimality_measures(problem, current_iterate);
    }
    
    /* generate the trial point */
    std::vector<double> trial_x = add_vectors(current_iterate.x, direction.x, step_length);
    Iterate trial_iterate(trial_x, direction.multipliers);
    double step_norm = step_length * direction.norm;
    this->subproblem.compute_optimality_measures(problem, trial_iterate);
    
    bool accept = false;
    /* check zero step */
    if (step_norm == 0.) {
        accept = true;
    }
    else {
        /* compute current exact l1 penalty: rho f + ||c|| */ 
        double current_exact_l1_penalty = direction.objective_multiplier * current_iterate.optimality_measure + current_iterate.feasibility_measure;
        double trial_exact_l1_penalty = direction.objective_multiplier * trial_iterate.optimality_measure + trial_iterate.feasibility_measure;
        
        /* check the validity of the trial step */
        double predicted_reduction = direction.predicted_reduction(problem, current_iterate, direction, step_length);
        double actual_reduction = current_exact_l1_penalty - trial_exact_l1_penalty;
        
        DEBUG << "Current: η = " << current_iterate.feasibility_measure << ", ω = " << current_iterate.optimality_measure << "\n";
        DEBUG << "Trial: η = " << trial_iterate.feasibility_measure << ", ω = " << trial_iterate.optimality_measure << "\n";
        DEBUG << "Predicted reduction: " << predicted_reduction << ", actual: " << actual_reduction << "\n\n";

        // Armijo sufficient decrease condition
        if (actual_reduction >= this->decrease_fraction_*predicted_reduction) {
            accept = true;
        }
    }
    
    if (accept) {
        statistics.add_statistic("penalty param.", direction.objective_multiplier);
        trial_iterate.compute_objective(problem);
        this->subproblem.compute_residuals(problem, trial_iterate, trial_iterate.multipliers, direction.objective_multiplier);
        DEBUG << "Residuals: ||c|| = " << trial_iterate.residuals.constraints << ", KKT = " << trial_iterate.residuals.KKT << ", cmpl = " << trial_iterate.residuals.complementarity << "\n";
        return std::optional<Iterate>{trial_iterate};
    }
    else {
        return std::nullopt;
    }
}