#include <cmath>
#include "TwoPhaseStrategy.hpp"
#include "Logger.hpp"

TwoPhaseStrategy::TwoPhaseStrategy(LocalApproximation& local_approximation, LocalSolutionConstants& constants, double tolerance):
		GlobalizationStrategy(local_approximation, tolerance), phase(OPTIMALITY), constants(constants) {
}

LocalSolution TwoPhaseStrategy::compute_step(Problem& problem, Iterate& current_point, double radius) {
	double objective_multiplier = (phase == OPTIMALITY) ? problem.obj_sign : 0.;
	LocalSolution solution = this->local_approximation.compute_optimality_step(problem, current_point, objective_multiplier, radius);
	solution.phase = OPTIMALITY;
	DEBUG << solution;

	if (solution.status == INFEASIBLE) { 
		/* infeasible subproblem during optimality phase */
		DEBUG << "Moving on to the feasibility step\n";
		
		/* different Hessian, but same constraint Jacobian */
		current_point.is_hessian_computed = false;
		/* partition of feasible and infeasible phase II constraints */
		ConstraintPartition constraint_partition = solution.constraint_partition;
		
		std::vector<double>& multipliers = (this->phase == OPTIMALITY) ? solution.multipliers : current_point.multipliers;
		/* compute the step in phase 1, starting from infeasible solution */
		solution = this->local_approximation.compute_infeasibility_step(problem, current_point, radius, solution.x, constraint_partition, multipliers);
		solution.phase = RESTORATION;
		solution.constraint_partition = constraint_partition;
		DEBUG << solution;
	}
	/* from this point on, the step is feasible */
	return solution;
}

void TwoPhaseStrategy::update_restoration_multipliers(Iterate& trial_point, ConstraintPartition& constraint_partition) {
	int number_variables = trial_point.x.size();
	for (unsigned int k = 0; k < constraint_partition.infeasible_set.size(); k++) {
		int j = constraint_partition.infeasible_set[k];
		if (constraint_partition.status[j] == INFEASIBLE_UPPER) {
			trial_point.multipliers[number_variables + j] = -1.;
		}
		else {
			trial_point.multipliers[number_variables + j] = 1.;
		}
	}
	return;
}

double TwoPhaseStrategy::compute_KKT_error(Problem& problem, Iterate& current_point) {
	double objective_multiplier = (this->phase == OPTIMALITY) ? 1. : 0.; 
	std::vector<double> lagrangian_gradient = this->compute_lagrangian_gradient(problem, current_point, objective_multiplier, current_point.multipliers);
	double KKTerror = norm_2(lagrangian_gradient);
	return KKTerror;
}
