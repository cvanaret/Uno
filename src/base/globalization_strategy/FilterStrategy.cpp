#include <iostream>
#include <cmath>
#include "FilterStrategy.hpp"
#include "Utils.hpp"

/* TODO option to switch between filter or non-monotonic filter */
FilterStrategy::FilterStrategy(LocalApproximation& local_approximation, Filter& filter_restoration,
	Filter& filter_optimality, LocalSolutionConstants& constants, Tolerances& tolerances, double tolerance):
		TwoPhaseStrategy(local_approximation, constants, tolerance), filter_restoration(filter_restoration), filter_optimality(filter_optimality), tolerances(tolerances) {
}

/* check acceptability of step(s) (filter & sufficient reduction)
 * precondition: feasible step
 * */
bool FilterStrategy::check_step(Problem& problem, Iterate& current_iterate, LocalSolution& solution, double step_length) {
	/* assemble trial point: TODO do not reevaluate if ||d|| = 0 */
	std::vector<double> x_trial = add_vectors(current_iterate.x, solution.x, step_length);
	Iterate trial_iterate(problem, x_trial, solution.multipliers);
	
	/* find out if transition of one phase to the other */
	if (this->phase == OPTIMALITY) {
		if (solution.phase == RESTORATION) {
			/* infeasible QP: go from phase II (optimality) to I (restoration) */
			DEBUG << "Switching from optimality to restoration phase\n";
			this->phase = RESTORATION;
			/* add [h,f] (c/s violation) to filter, entering restoration */
			this->filter_optimality.add(current_iterate.infeasibility_measure, current_iterate.feasibility_measure);
			
			/* re-initialize the restoration filter */
			current_iterate.infeasibility_measure = problem.feasible_residual_norm(solution.constraint_partition, current_iterate.constraints);
			current_iterate.feasibility_measure = problem.infeasible_residual_norm(solution.constraint_partition, current_iterate.constraints);
			current_iterate.is_hessian_computed = false;
			this->filter_restoration.reset();
			this->filter_restoration.upper_bound = this->filter_optimality.upper_bound;
			this->filter_restoration.add(current_iterate.infeasibility_measure, current_iterate.feasibility_measure);
		}
	}
	else { /* phase == RESTORATION */
		/* check whether we can switch from phase I (restoration) to II (optimality) */
		if (solution.phase == OPTIMALITY && this->filter_optimality.query(trial_iterate.infeasibility_measure, trial_iterate.feasibility_measure)) {
			DEBUG << "Switching from restoration to optimality phase\n";
			this->phase = OPTIMALITY;
			current_iterate.infeasibility_measure = current_iterate.residual;
			current_iterate.feasibility_measure = current_iterate.objective;
		}
	}
	
	/* if RESTORATION phase, compute (in)feasibility measures of trial point */
	if (this->phase == RESTORATION) {
		trial_iterate.infeasibility_measure = problem.feasible_residual_norm(solution.constraint_partition, trial_iterate.constraints);
		trial_iterate.feasibility_measure = problem.infeasible_residual_norm(solution.constraint_partition, trial_iterate.constraints);
	}
	
	DEBUG << "TESTING " << current_iterate.infeasibility_measure << " " << current_iterate.feasibility_measure << " ";
	DEBUG << trial_iterate.infeasibility_measure << " " << trial_iterate.feasibility_measure << "\n";
	
	Filter& filter = (this->phase == OPTIMALITY) ? this->filter_optimality : this->filter_restoration;
	
	/* check acceptance */
	bool accept = false;
	bool acceptable = filter.query(trial_iterate.infeasibility_measure, trial_iterate.feasibility_measure);
	if (acceptable) {
		// check acceptance wrt current x (h,f)
		acceptable = filter.query_current_iterate(current_iterate.infeasibility_measure, current_iterate.feasibility_measure, trial_iterate.infeasibility_measure, trial_iterate.feasibility_measure);
		if (acceptable) {
			double predicted_reduction = -solution.objective;
			double actual_reduction = filter.compute_actual_reduction(current_iterate.feasibility_measure, current_iterate.infeasibility_measure, trial_iterate.feasibility_measure);
			DEBUG << "Predicted reduction: " << predicted_reduction << ", actual: " << actual_reduction << "\n\n";
			
			if (predicted_reduction < this->constants.Delta*current_iterate.infeasibility_measure*current_iterate.infeasibility_measure) {
				filter.add(current_iterate.infeasibility_measure, current_iterate.feasibility_measure);
				accept = true;
			}
			else if (actual_reduction >= this->constants.Sigma*std::max(0., predicted_reduction-1e-9)) {
				accept = true;
			}
		}
	}
	
	/* correct multipliers for infeasibility problem */
	if (accept) {
		if(solution.phase == RESTORATION) {
			this->update_restoration_multipliers(trial_iterate, solution.constraint_partition);
		}
		current_iterate = trial_iterate;
		current_iterate.KKTerror = this->compute_KKT_error(problem, current_iterate);
		current_iterate.complementarity_error = (this->phase == OPTIMALITY) ? this->compute_complementarity_error(problem, current_iterate) : 0.;
		double step_norm = step_length*norm_inf(solution.x);
		current_iterate.status = this->compute_status(problem, current_iterate, step_norm);
		INFO << "phase: " << this->phase << "\t";
	}
	
	return accept;
}

OptimalityStatus FilterStrategy::compute_status(Problem& problem, Iterate& current_iterate, double step_norm) {
	OptimalityStatus status = NOT_OPTIMAL;
	
	if (current_iterate.residual <= this->tolerance*problem.number_constraints) {
		if (current_iterate.KKTerror <= this->tolerance*sqrt(problem.number_variables) &&
				current_iterate.complementarity_error <= this->tolerance*(problem.number_variables + problem.number_constraints)) {
			if (this->phase == OPTIMALITY) {
				status = KKT_POINT;
			}
			else {
				status = FJ_POINT;
			}
		}
		else if (step_norm <= this->tolerance/100.) {
			status = FEASIBLE_SMALL_STEP;
		}
	}
	else if (step_norm <= this->tolerance/100.) {
		status = INFEASIBLE_SMALL_STEP;
	}
	return status;
}

void FilterStrategy::initialize(Problem& problem, Iterate& current_iterate) {
	/* set the filter upper bound */
	double upper_bound = std::max(this->tolerances.ubd, this->tolerances.fact*current_iterate.residual);
	this->filter_optimality.upper_bound = upper_bound;
	this->filter_restoration.upper_bound = upper_bound;
	
	/* allocate the subproblem solver */
	this->local_approximation.allocate_solver(problem.number_variables, problem.number_constraints);
	return;
}
