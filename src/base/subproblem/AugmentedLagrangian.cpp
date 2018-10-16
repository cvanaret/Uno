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
    
    /* identify the inequality constraint slacks */
    std::map<int,int> slacked_constraints;
    int current_slack = 0;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] != EQUAL_BOUNDS) {
            slacked_constraints[j] = current_slack;
            current_slack++;
            // add slack as a primal variable
            first_iterate.x.push_back(first_iterate.constraints[j]);
        }
    }
    
    std::cout << "SLACKED CONSTRAINTS:\n";
    for (std::pair<const int, int> element: slacked_constraints) {
        std::cout << "c" << element.first << " has slack s" << element.second << std::endl;  // ? or is abc an iterator?
    }
    
    /* TODO allocate the solver */
    this->solver.initialize(slacked_constraints);
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
