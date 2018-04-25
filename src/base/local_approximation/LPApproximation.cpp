#include <iostream>
#include <map>
#include "LPApproximation.hpp"
#include "Vector.hpp"

LPApproximation::LPApproximation(LPSolver& solver): LocalApproximation("LP"), solver(solver) {
}

LocalSolution LPApproximation::compute_direction(Problem& problem, Phase& phase, Iterate& current_point) {
	// TODO check phase and generate the right suproblem
	LP lp = this->generate_LP(problem, current_point);
	std::vector<double> d0(lp.number_variables); // = {0.}
	LocalSolution solution = this->solver.solve(lp, d0);
	return solution;
}

LP LPApproximation::generate_LP(Problem& problem, Iterate& current_point) const {
	LP lp(problem.number_variables, problem.number_constraints);

	std::map<int,double> objective_jacobian = problem.objective_sparse_gradient(current_point.x);
	current_point.set_objective_jacobian(objective_jacobian);
	/* Jacobian and its sparsity pattern */
	//lp.jacobian = current_point.compute_jacobian(problem);
	lp.jacobian_sparsity = problem.jacobian_sparsity;

	/* variables range intersected with trust region */
	for (int i = 0; i < lp.number_variables; i++) {
		lp.lb[i] = problem.variables[i].lb - current_point.x[i];
		lp.ub[i] = problem.variables[i].ub - current_point.x[i];
	}
	/* linearized constraints */
	for (int j = 0; j < lp.number_constraints; j++) {
		lp.lb[problem.number_variables + j] = problem.constraints[j].lb - current_point.constraints[j];
		lp.ub[problem.number_variables + j] = problem.constraints[j].ub - current_point.constraints[j];
	}

	return lp;
}
