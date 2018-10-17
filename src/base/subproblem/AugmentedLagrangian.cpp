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
    
    /* initialize the solver */
    this->solver.initialize(slacked_constraints);
    //this->solver.allocate(number_variables, number_constraints);
    return first_iterate;
}

LocalSolution AugmentedLagrangian::compute_optimality_step(Problem& problem, Iterate& current_iterate, double radius) {
    LocalSolution solution = this->solver.solve(problem, current_iterate);
    /* the Augmented Lagrangian subproblem returns a new iterate. We should now compute the step */
    std::vector<double> step(solution.x);
    for (unsigned int i = 0; i < current_iterate.x.size(); i++) {
        step[i] -= current_iterate.x[i];
    }
    solution.x = step;
    solution.phase_1_required = this->phase_1_required(solution);
    return solution;
}

LocalSolution AugmentedLagrangian::compute_infeasibility_step(Problem& problem, Iterate& current_iterate, double radius, LocalSolution& phase_II_solution) {
    std::cout << "AugmentedLagrangian::compute_infeasibility_step not implemented yet\n";
    throw std::out_of_range("Not implemented yet.");
}

LocalSolution AugmentedLagrangian::compute_l1_penalty_step(Problem& problem, Iterate& current_iterate, double radius, double penalty_parameter, PenaltyDimensions penalty_dimensions) {
    std::cout << "AugmentedLagrangian::compute_l1_penalty_step not implemented yet\n";
    throw std::out_of_range("Not implemented yet.");
}

void AugmentedLagrangian::compute_measures(Problem& problem, Iterate& iterate) {
    return;
}

bool AugmentedLagrangian::phase_1_required(LocalSolution& solution) {
    return solution.status == INFEASIBLE;
}