#include <iostream>
#include <map>
#include "LPApproximation.hpp"
#include "Vector.hpp"

LPApproximation::LPApproximation(LPSolver& solver): LocalApproximation("LP"), solver(solver) {
}

SubproblemSolution LPApproximation::compute_direction(Problem& problem, Phase& phase, Iterate& current_iterate) {
	// TODO check phase and generate the right suproblem
	LP lp = this->generate_LP(problem, current_iterate);
	std::vector<double> d0(lp.number_variables); // = {0.}
	SubproblemSolution solution = this->solver.solve(lp, d0);
	return solution;
}

LP LPApproximation::generate_LP(Problem& problem, Iterate& current_iterate) const {
	LP lp(problem.number_variables, problem.number_constraints);

	std::map<int,double> objective_jacobian = problem.objective_sparse_gradient(current_iterate.x);
	current_iterate.set_objective_jacobian(objective_jacobian);
	/* Jacobian and its sparsity pattern */
	//lp.jacobian = current_iterate.compute_jacobian(problem);
	lp.jacobian_sparsity = problem.jacobian_sparsity;

	/* variables range intersected with trust region */
	for (int i = 0; i < lp.number_variables; i++) {
		lp.lb[i] = problem.variables[i].lb - current_iterate.x[i];
		lp.ub[i] = problem.variables[i].ub - current_iterate.x[i];
	}
	/* linearized constraints */
	for (int j = 0; j < lp.number_constraints; j++) {
		lp.lb[problem.number_variables + j] = problem.constraints[j].lb - current_iterate.constraints[j];
		lp.ub[problem.number_variables + j] = problem.constraints[j].ub - current_iterate.constraints[j];
	}

	return lp;
}
