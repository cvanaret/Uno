#include <iostream>
#include <cmath>
#include "Argonot.hpp"
#include "FilterStrategy.hpp"
#include "Utils.hpp"

/* TODO option to switch between filter or non-monotonic filter */
FilterStrategy::FilterStrategy(Subproblem& subproblem, std::shared_ptr<Filter> filter_optimality, std::shared_ptr<Filter> filter_restoration, FilterStrategyParameters& strategy_parameters, double tolerance):
GlobalizationStrategy(subproblem, tolerance), filter_optimality(filter_optimality), filter_restoration(filter_restoration), current_phase(OPTIMALITY), parameters(strategy_parameters) {
}

Iterate FilterStrategy::initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers) {
    /* initialize the subproblem */
    Iterate first_iterate = this->subproblem.evaluate_initial_point(problem, x, multipliers);
    this->subproblem.compute_residuals(problem, first_iterate, first_iterate.multipliers, 1.);

    /* set the filter upper bound */
    double upper_bound = std::max(this->parameters.ubd, this->parameters.fact * first_iterate.feasibility_measure);
    this->filter_optimality->upper_bound = upper_bound;
    this->filter_restoration->upper_bound = upper_bound;
    return first_iterate;
}

/* check acceptability of step(s) (filter & sufficient reduction)
 * precondition: feasible step
 * */
bool FilterStrategy::check_step(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution, double step_length) {
    /* check if subproblem definition changed */
    if (this->subproblem.subproblem_definition_changed) {
        this->filter_optimality->reset();
        this->subproblem.subproblem_definition_changed = false;
        this->subproblem.compute_optimality_measures(problem, current_iterate);
    }

    /* assemble trial point: TODO do not reevaluate if ||d|| = 0 */
    std::vector<double> trial_x = add_vectors(current_iterate.x, solution.x, step_length);
    Iterate trial_iterate(trial_x, solution.multipliers);
    double step_norm = step_length * solution.norm;
    this->subproblem.compute_optimality_measures(problem, trial_iterate);
    
    bool accept = false;
    /* check zero step */
    if (step_norm == 0.) {
        accept = true;
    }
    else {
        /* possibly switch phase */
        this->switch_phase(problem, solution, current_iterate, trial_iterate);

        /* if RESTORATION phase, compute (in)feasibility measures of trial point */
        if (this->current_phase == RESTORATION) {
            this->subproblem.compute_infeasibility_measures(problem, trial_iterate, solution);
        }

        DEBUG << "Current: η = " << current_iterate.feasibility_measure << ", ω = " << current_iterate.optimality_measure << "\n";
        DEBUG << "Trial: η = " << trial_iterate.feasibility_measure << ", ω = " << trial_iterate.optimality_measure << "\n";

        /* check acceptance */
        Filter& filter = (this->current_phase == OPTIMALITY) ? *(this->filter_optimality) : *(this->filter_restoration);
        bool acceptable = filter.accept(trial_iterate.feasibility_measure, trial_iterate.optimality_measure);
        if (acceptable) {
            // check acceptance wrt current x (h,f)
            acceptable = filter.improves_current_iterate(current_iterate.feasibility_measure, current_iterate.optimality_measure, trial_iterate.feasibility_measure, trial_iterate.optimality_measure);
            if (acceptable) {
                double predicted_reduction = this->subproblem.compute_predicted_reduction(problem, current_iterate, solution, step_length);
                double actual_reduction = filter.compute_actual_reduction(current_iterate.optimality_measure, current_iterate.feasibility_measure, trial_iterate.optimality_measure);
                DEBUG << "Predicted reduction: " << predicted_reduction << ", actual: " << actual_reduction << "\n\n";

                /* switching condition */
                if (predicted_reduction < this->parameters.Delta * std::pow(current_iterate.feasibility_measure, 2)) {
                    filter.add(current_iterate.feasibility_measure, current_iterate.optimality_measure);
                    accept = true;
                }
                /* Armijo sufficient decrease condition: predicted_reduction should be positive */
                else if (actual_reduction >= this->parameters.Sigma*step_length*std::max(0., predicted_reduction - 1e-9)) {
                //else if (predicted_reduction > 0. && actual_reduction >= this->parameters.Sigma*predicted_reduction) {
                    accept = true;
                }
            }
        }
    }

    /* correct multipliers for infeasibility problem */
    if (accept) {
        if (solution.phase == RESTORATION) {
            this->update_restoration_multipliers(trial_iterate, solution.constraint_partition);
        }
        trial_iterate.compute_objective(problem);
        this->subproblem.compute_residuals(problem, trial_iterate, trial_iterate.multipliers, solution.objective_multiplier);
        trial_iterate.status = this->compute_status(problem, trial_iterate, step_norm, solution.objective_multiplier);
        INFO << "phase: " << this->current_phase << "\t";
        current_iterate = trial_iterate;
        DEBUG << "Residuals: ||c|| = " << current_iterate.residuals.constraints << ", KKT = " << current_iterate.residuals.KKT << ", cmpl = " << current_iterate.residuals.complementarity << "\n";
    }
    return accept;
}

void FilterStrategy::switch_phase(Problem& problem, SubproblemSolution& solution, Iterate& current_iterate, Iterate& trial_iterate) {
    /* find out if transition of one phase to the other */
    if (this->current_phase == OPTIMALITY) {
        if (solution.phase == RESTORATION) {
            /* infeasible QP: go from phase II (optimality) to I (restoration) */
            DEBUG << "Switching from optimality to restoration phase\n";
            this->current_phase = RESTORATION;
            /* add [h,f] (c/s violation) to filter, entering restoration */
            this->filter_optimality->add(current_iterate.feasibility_measure, current_iterate.optimality_measure);

            /* re-initialize the restoration filter */
            this->filter_restoration->reset();
            this->filter_restoration->upper_bound = this->filter_optimality->upper_bound;
            this->subproblem.compute_infeasibility_measures(problem, current_iterate, solution);
            this->filter_restoration->add(current_iterate.feasibility_measure, current_iterate.optimality_measure);
            current_iterate.is_hessian_computed = false;
        }
    }
        /* check whether we can switch from phase I (restoration) to II (optimality) */
    else if (solution.phase == OPTIMALITY && this->filter_optimality->accept(trial_iterate.feasibility_measure, trial_iterate.optimality_measure)) {
        DEBUG << "Switching from restoration to optimality phase\n";
        this->current_phase = OPTIMALITY;
        this->subproblem.compute_optimality_measures(problem, current_iterate);
    }
    return;
}

void FilterStrategy::update_restoration_multipliers(Iterate& trial_iterate, ConstraintPartition& constraint_partition) {
    for (int j: constraint_partition.infeasible) {
        if (constraint_partition.constraint_feasibility[j] == INFEASIBLE_UPPER) {
            trial_iterate.multipliers.constraints[j] = -1.;
        }
        else {
            trial_iterate.multipliers.constraints[j] = 1.;
        }
    }
    return;
}