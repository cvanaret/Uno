#include <iostream>
#include <cmath>
#include "Argonot.hpp"
#include "FilterStrategy.hpp"
#include "Utils.hpp"

/* TODO option to switch between filter or non-monotonic filter */
FilterStrategy::FilterStrategy(Subproblem& subproblem, std::shared_ptr<Filter> filter_optimality, std::shared_ptr<Filter> filter_restoration, FilterStrategyConstants& constants, double tolerance):
GlobalizationStrategy(subproblem, tolerance), filter_optimality(filter_optimality), filter_restoration(filter_restoration), current_phase(OPTIMALITY), constants(constants) {
}

Iterate FilterStrategy::initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers, bool use_trust_region) {
    /* initialize the subproblem */
    Iterate first_iterate = this->subproblem.initialize(problem, x, multipliers, use_trust_region);

    first_iterate.KKT_residual = Argonot::compute_KKT_error(problem, first_iterate, 1.);
    first_iterate.complementarity_residual = Argonot::compute_complementarity_error(problem, first_iterate);

    /* set the filter upper bound */
    double upper_bound = std::max(this->constants.ubd, this->constants.fact * first_iterate.feasibility_measure);
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
    this->subproblem.compute_optimality_measures(problem, trial_iterate);
    double step_norm = step_length * solution.norm;

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

        DEBUG << "TESTING trial (" << trial_iterate.feasibility_measure << ", " << trial_iterate.optimality_measure << ") against current (";
        DEBUG << current_iterate.feasibility_measure << ", " << current_iterate.optimality_measure << ")\n";

        /* check acceptance */
        Filter& filter = (this->current_phase == OPTIMALITY) ? *(this->filter_optimality) : *(this->filter_restoration);
        bool acceptable = filter.accept(trial_iterate.feasibility_measure, trial_iterate.optimality_measure);
        if (acceptable) {
            // check acceptance wrt current x (h,f)
            acceptable = filter.improves_current_iterate(current_iterate.feasibility_measure, current_iterate.optimality_measure, trial_iterate.feasibility_measure, trial_iterate.optimality_measure);
            if (acceptable) {
                double predicted_reduction = this->subproblem.compute_predicted_reduction(current_iterate, solution);
                double actual_reduction = filter.compute_actual_reduction(current_iterate.optimality_measure, current_iterate.feasibility_measure, trial_iterate.optimality_measure);
                DEBUG << "Predicted reduction: " << step_length*predicted_reduction << ", actual: " << actual_reduction << "\n\n";

                /* switching condition */
                if (step_length*predicted_reduction < this->constants.Delta * std::pow(current_iterate.feasibility_measure, 2)) {
                    filter.add(current_iterate.feasibility_measure, current_iterate.optimality_measure);
                    accept = true;
                }
                    /* Armijo sufficient decrease condition */
                // else if (actual_reduction >= this->constants.Sigma*step_length*std::max(0., predicted_reduction - 1e-9)) {
                else if (actual_reduction >= this->constants.Sigma*step_length*predicted_reduction) {
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
        trial_iterate.compute_constraint_residual(problem, this->subproblem.residual_norm);
        trial_iterate.KKT_residual = Argonot::compute_KKT_error(problem, trial_iterate, solution.objective_multiplier);
        trial_iterate.complementarity_residual = (this->current_phase == OPTIMALITY) ? Argonot::compute_complementarity_error(problem, trial_iterate) : 0.;
        trial_iterate.status = this->compute_status(problem, trial_iterate, step_norm, solution.objective_multiplier);
        INFO << "phase: " << this->current_phase << "\t";
        current_iterate = trial_iterate;
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

OptimalityStatus FilterStrategy::compute_status(Problem& problem, Iterate& current_iterate, double step_norm, double objective_multiplier) {
    OptimalityStatus status = NOT_OPTIMAL;

    /* TODO: check if test on residual can indeed be replaced by infeasibility_measure */
    if (current_iterate.constraint_residual <= this->tolerance * problem.number_constraints) {
        if (current_iterate.KKT_residual <= this->tolerance * std::sqrt(problem.number_variables) && current_iterate.complementarity_residual <= this->tolerance * (problem.number_variables + problem.number_constraints)) {
            if (this->current_phase == OPTIMALITY) {
                status = KKT_POINT;
            }
            else {
                status = FJ_POINT;
            }
        }
        else if (step_norm <= this->tolerance / 100.) {
            status = FEASIBLE_SMALL_STEP;
        }
    }
    else if (step_norm <= this->tolerance / 100.) {
        status = INFEASIBLE_SMALL_STEP;
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