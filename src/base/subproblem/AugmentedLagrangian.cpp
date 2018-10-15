#include <iostream>
#include <cmath>
#include "AugmentedLagrangian.hpp"
#include "Utils.hpp"

AugmentedLagrangian::AugmentedLagrangian(): Subproblem() {
}

Iterate AugmentedLagrangian::initialize(Problem& problem, std::vector<double>& x, std::vector<double>& bound_multipliers, std::vector<double>& constraint_multipliers, int number_variables, int number_constraints, bool use_trust_region) {
    Iterate first_iterate(problem, x, bound_multipliers, constraint_multipliers);

    /* compute the optimality and feasibility measures of the initial point */
    this->compute_measures(problem, first_iterate);

    /* allocate the solver */
    //this->solver.allocate(number_variables, number_constraints);
    return first_iterate;
}

LocalSolution AugmentedLagrangian::compute_optimality_step(Problem& problem, Iterate& current_iterate, double radius) {
	std::vector<double> x;
	LocalSolution solution = this->solver.solve(problem, current_iterate);
	return solution;
}

LocalSolution AugmentedLagrangian::compute_infeasibility_step(Problem& problem, Iterate& current_iterate, double radius, LocalSolution& phase_II_solution) {
	std::vector<double> x;
	LocalSolution solution(x, x, x);
	return solution;
}

LocalSolution AugmentedLagrangian::compute_l1_penalty_step(Problem& problem, Iterate& current_iterate, double radius, double penalty_parameter, PenaltyDimensions penalty_dimensions) {
	std::vector<double> x;
	LocalSolution solution(x, x, x);
	return solution;
}

void AugmentedLagrangian::compute_measures(Problem& problem, Iterate& iterate) {
	return;
}
