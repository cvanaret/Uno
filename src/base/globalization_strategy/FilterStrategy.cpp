#include <iostream>
#include <cmath>
#include "FilterStrategy.hpp"
#include "Utils.hpp"

/* TODO option to switch between filter or non-monotonic filter */
FilterStrategy::FilterStrategy(Subproblem& subproblem, std::shared_ptr<Filter> filter_optimality,
                               std::shared_ptr<Filter> filter_restoration, TwoPhaseConstants& constants, double tolerance) :
TwoPhaseStrategy(subproblem, constants, tolerance), filter_optimality(filter_optimality), filter_restoration(filter_restoration) {
}

Iterate FilterStrategy::initialize(Problem& problem, std::vector<double>& x, std::vector<double>& bound_multipliers, std::vector<double>& constraint_multipliers, bool use_trust_region) {
    /* initialize the subproblem */
    Iterate first_iterate = this->subproblem.initialize(problem, x, bound_multipliers, constraint_multipliers, problem.number_variables, problem.number_constraints, use_trust_region);

    first_iterate.KKTerror = this->compute_KKT_error(problem, first_iterate, 1.);
    first_iterate.complementarity_error = this->compute_complementarity_error(problem, first_iterate);

    /* set the filter upper bound */
    double upper_bound = std::max(this->constants.ubd, this->constants.fact * first_iterate.feasibility_measure);
    this->filter_optimality->upper_bound = upper_bound;
    this->filter_restoration->upper_bound = upper_bound;
    return first_iterate;
}

/* check acceptability of step(s) (filter & sufficient reduction)
 * precondition: feasible step
 * */
bool FilterStrategy::check_step(Problem& problem, Iterate& current_iterate, LocalSolution& solution, double step_length) {
    /* assemble trial point: TODO do not reevaluate if ||d|| = 0 */
    std::vector<double> trial_x = add_vectors(current_iterate.x, solution.x, step_length);
    Iterate trial_iterate(problem, trial_x, solution.bound_multipliers, solution.constraint_multipliers);
    this->compute_measures(problem, trial_iterate);
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
        if (this->phase == RESTORATION) {
            trial_iterate.feasibility_measure = problem.feasible_residual_norm(solution.constraint_partition, trial_iterate.constraints);
            trial_iterate.optimality_measure = problem.infeasible_residual_norm(solution.constraint_partition, trial_iterate.constraints);
        }

        DEBUG << "TESTING trial " << trial_iterate.feasibility_measure << " " << trial_iterate.optimality_measure << " against current ";
        DEBUG << current_iterate.feasibility_measure << " " << current_iterate.optimality_measure << "\n";

        /* check acceptance */
        Filter& filter = (this->phase == OPTIMALITY) ? *(this->filter_optimality) : *(this->filter_restoration);
        bool acceptable = filter.query(trial_iterate.feasibility_measure, trial_iterate.optimality_measure);
        if (acceptable) {
            // check acceptance wrt current x (h,f)
            acceptable = filter.query_current_iterate(current_iterate.feasibility_measure, current_iterate.optimality_measure, trial_iterate.feasibility_measure, trial_iterate.optimality_measure);
            if (acceptable) {
                double predicted_reduction = this->compute_predicted_reduction(solution, step_length);
                double actual_reduction = filter.compute_actual_reduction(current_iterate.optimality_measure, current_iterate.feasibility_measure, trial_iterate.optimality_measure);
                DEBUG << "Predicted reduction: " << predicted_reduction << ", actual: " << current_iterate.optimality_measure << "-" << trial_iterate.optimality_measure << " = " << actual_reduction << "\n\n";

                /* switching condition */
                if (predicted_reduction < this->constants.Delta * std::pow(current_iterate.feasibility_measure, 2)) {
                    filter.add(current_iterate.feasibility_measure, current_iterate.optimality_measure);
                    accept = true;
                }
                    /* sufficient decrease condition */
                    //else if (actual_reduction >= this->constants.Sigma * std::max(0., predicted_reduction - 1e-9)) {
                else if (actual_reduction >= this->constants.Sigma * predicted_reduction) {
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
        current_iterate = trial_iterate;
        double objective_multiplier = (this->phase == OPTIMALITY) ? 1. : 0.;
        current_iterate.KKTerror = this->compute_KKT_error(problem, current_iterate, objective_multiplier);
        current_iterate.complementarity_error = (this->phase == OPTIMALITY) ? this->compute_complementarity_error(problem, current_iterate) : 0.;
        current_iterate.status = this->compute_status(problem, current_iterate, step_norm);
        INFO << "phase: " << this->phase << "\t";
    }

    return accept;
}

/* compute the predicted reduction, taking into account the step length */
double FilterStrategy::compute_predicted_reduction(LocalSolution& solution, double step_length) {
    if (step_length == 1.) {
        return -solution.objective;
    }
    else {
        return -(step_length * solution.objective_terms.linear + step_length * step_length * solution.objective_terms.quadratic);
    }
}

void FilterStrategy::switch_phase(Problem& problem, LocalSolution& solution, Iterate& current_iterate, Iterate& trial_iterate) {
    /* find out if transition of one phase to the other */
    if (this->phase == OPTIMALITY) {
        if (solution.phase == RESTORATION) {
            /* infeasible QP: go from phase II (optimality) to I (restoration) */
            DEBUG << "Switching from optimality to restoration phase\n";
            this->phase = RESTORATION;
            /* add [h,f] (c/s violation) to filter, entering restoration */
            this->filter_optimality->add(current_iterate.feasibility_measure, current_iterate.optimality_measure);

            /* re-initialize the restoration filter */
            current_iterate.feasibility_measure = problem.feasible_residual_norm(solution.constraint_partition, current_iterate.constraints);
            current_iterate.optimality_measure = problem.infeasible_residual_norm(solution.constraint_partition, current_iterate.constraints);
            current_iterate.is_hessian_computed = false;
            this->filter_restoration->reset();
            this->filter_restoration->upper_bound = this->filter_optimality->upper_bound;
            this->filter_restoration->add(current_iterate.feasibility_measure, current_iterate.optimality_measure);
        }
    }
        /* check whether we can switch from phase I (restoration) to II (optimality) */
    else if (solution.phase == OPTIMALITY && this->filter_optimality->query(trial_iterate.feasibility_measure, trial_iterate.optimality_measure)) {
        DEBUG << "Switching from restoration to optimality phase\n";
        this->phase = OPTIMALITY;
        this->compute_measures(problem, current_iterate);
    }
    return;
}

OptimalityStatus FilterStrategy::compute_status(Problem& problem, Iterate& current_iterate, double step_norm) {
    OptimalityStatus status = NOT_OPTIMAL;

    /* TODO: check if test on residual can indeed be replaced by infeasibility_measure */
    if (current_iterate.residual <= this->tolerance * problem.number_constraints) {
        if (current_iterate.KKTerror <= this->tolerance * std::sqrt(problem.number_variables) &&
                current_iterate.complementarity_error <= this->tolerance * (problem.number_variables + problem.number_constraints)) {
            if (this->phase == OPTIMALITY) {
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
    return status;
}

void FilterStrategy::compute_measures(Problem& problem, Iterate& iterate) {
    this->subproblem.compute_measures(problem, iterate);
    return;
}
