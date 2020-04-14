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
    std::map<int,int> slacked_constraints;                     // sparse dictionary of inequality c/s
    int current_slack = 0;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] != EQUAL_BOUNDS) {
            slacked_constraints[j] = current_slack;
            current_slack++;
            // add slack as a primal variable
            first_iterate.x.push_back(first_iterate.constraints[j]);
            //first_iterate.x.push_back(problem.constraint_lb[j]);
        }
    }
    
    /* initialize the solver */
    this->solver.initialize(slacked_constraints);
    //this->solver.allocate(number_variables, number_constraints);
    return first_iterate;
}

//// evaluate the augmented Lagrangian at x,y, where y=constraint_multipliers
//double compute_augmented_lagrangian_(Problem& problem, std::vector<double>& x, std::vector<double>& constraints, std::vector<double>& constraint_multipliers) {
    //// contribution of the objective
    //double f = problem.objective(x);
    //// contribution of the constraints
    //for (int j = 0; j < problem.number_constraints; j++) {
        //double constraint_value;
        //try {
            //// inequality constraint: need to subtract slack values
            //int current_slack = this->slacked_constraints_.at(j);
            //constraint_value = constraints[j] - x[problem.number_variables + current_slack];
        //}
        //catch (std::out_of_range) {
            //// equality constraint
            //constraint_value = constraints[j] - problem.constraint_lb[j];
        //}
        //f -= constraint_multipliers[j]*constraint_value;       // f = f - lambda[i]*c[i]
        //f += this->rho/2.*constraint_value*constraint_value;   // f = f + rho/2*(c[i])^2 = augmented Lagrangian
    //}
    //return f;
//}

//// evaluate the gradient of the augmented Lagrangian at x,y, where y=constraint_multipliers
//std::vector<double> compute_augmented_lagrangian_gradient_(Problem& problem, std::vector<double>& x, std::vector<double>& constraints, std::vector<double>& constraint_multipliers) {
    //// start with gradient of the objective
    //std::vector<double> augmented_lagrangian_gradient = problem.objective_dense_gradient(x);
    //// gradient of the constraints wrt the variables
    //for (int j = 0; j < problem.number_constraints; j++) {
        //double constraint_value;
        //try {
            //// inequality constraint: need to subtract slacks
            //int current_slack = this->slacked_constraints_.at(j);
            //constraint_value = constraints[j] - x[problem.number_variables + current_slack];
        //}
        //catch (std::out_of_range) {
            //// equality constraint
            //constraint_value = constraints[j] - problem.constraint_lb[j];
        //}
        //double factor = constraint_multipliers[j] - this->rho*constraint_value;
        //// add the gradient contribution from the constraints
        //std::vector<double> constraint_gradient = problem.constraint_dense_gradient(j, x);
        //for (int i = 0; i < problem.number_variables; i++) {
            //augmented_lagrangian_gradient[i] -= factor*constraint_gradient[i];
        //}
    //}
    //// gradient of the constraints wrt the slacks
    //for (std::pair<const int, int> element: slacked_constraints_) {
        //int j = element.first;                // index of the constraint
        //int current_slack = element.second;   // index of the slack in [0, number_of_slacks[
        //double constraint_value = constraints[j] - x[problem.number_variables + current_slack];
        //double derivative = constraint_multipliers[j] - this->rho*constraint_value;
        //augmented_lagrangian_gradient.push_back(derivative); // sticks gradient terms at end of n (number of vars) grad.
    //}
    //return augmented_lagrangian_gradient;
//}

SubproblemSolution AugmentedLagrangian::compute_optimality_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds) {
    /* Solve (approx/exact) the augm. Lagrangian subproblem */ 
    SubproblemSolution solution = this->solver.solve(problem, current_iterate);
    /* the Augmented Lagrangian subproblem returns a new iterate. We should now compute the step */
    std::vector<double> step(solution.x); /* copy/contruct step := solution.x */
    for (unsigned int i = 0; i < current_iterate.x.size(); i++) {
        step[i] -= current_iterate.x[i];
    }
    solution.x = step;
    solution.phase_1_required = this->phase_1_required(solution);
    return solution;
}

SubproblemSolution AugmentedLagrangian::compute_infeasibility_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds, SubproblemSolution& phase_II_solution) {
    std::cout << "AugmentedLagrangian::compute_infeasibility_step not implemented yet\n";
    throw std::out_of_range("Not implemented yet.");
}

//SubproblemSolution AugmentedLagrangian::compute_l1_penalty_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds, double penalty_parameter, PenaltyDimensions penalty_dimensions) {
//    std::cout << "AugmentedLagrangian::compute_l1_penalty_step not implemented yet\n";
//    throw std::out_of_range("Not implemented yet.");
//}

void AugmentedLagrangian::compute_measures(Problem& problem, Iterate& iterate) {
    return;
}

bool AugmentedLagrangian::phase_1_required(SubproblemSolution& solution) {
    return solution.status == INFEASIBLE;
}
