#include <cmath>
#include "TwoPhaseStrategy.hpp"
#include "Logger.hpp"

TwoPhaseStrategy::TwoPhaseStrategy(LocalApproximation& local_approximation, TwoPhaseConstants& constants, double tolerance):
		GlobalizationStrategy(local_approximation, tolerance), phase(OPTIMALITY), constants(constants) {
}

LocalSolution TwoPhaseStrategy::compute_step(Problem& problem, Iterate& current_iterate, double radius) {
	double objective_multiplier = problem.obj_sign;
	//double objective_multiplier = (phase == OPTIMALITY) ? problem.obj_sign : 0.;
	LocalSolution solution = this->local_approximation.compute_optimality_step(problem, current_iterate, objective_multiplier, radius);
	solution.phase = OPTIMALITY;
	DEBUG << solution;
	
	if (solution.status == INFEASIBLE) { 
		/* infeasible subproblem during optimality phase */
		DEBUG << "Moving on to the feasibility step\n";
		
		/* different Hessian, but same constraint Jacobian */
		current_iterate.is_hessian_computed = false;
		/* partition of feasible and infeasible phase II constraints */
		ConstraintPartition constraint_partition = solution.constraint_partition;
		
		/* compute the step in phase 1, starting from infeasible solution */
		solution = this->local_approximation.compute_infeasibility_step(problem, current_iterate, radius, solution.x, constraint_partition, solution.multipliers);
		solution.constraint_partition = constraint_partition;
		solution.phase = RESTORATION;
		DEBUG << solution;
	}
	/* from this point on, the step is feasible */
	return solution;
}

void TwoPhaseStrategy::update_restoration_multipliers(Iterate& trial_iterate, ConstraintPartition& constraint_partition) {
	int number_variables = trial_iterate.x.size();
	for (unsigned int k = 0; k < constraint_partition.infeasible_set.size(); k++) {
		int j = constraint_partition.infeasible_set[k];
		if (constraint_partition.status[j] == INFEASIBLE_UPPER) {
			trial_iterate.multipliers[number_variables + j] = -1.;
		}
		else {
			trial_iterate.multipliers[number_variables + j] = 1.;
		}
	}
	return;
}

double TwoPhaseStrategy::compute_KKT_error(Problem& problem, Iterate& current_iterate) {
	double objective_multiplier = (this->phase == OPTIMALITY) ? 1. : 0.; 
	std::vector<double> lagrangian_gradient = this->compute_lagrangian_gradient(problem, current_iterate, objective_multiplier, current_iterate.multipliers);
	double KKTerror = norm_2(lagrangian_gradient);
	return KKTerror;
}
