#include <iostream>
#include <vector>
#include "Logger.hpp"
#include "Filter.hpp"
#include "LBFGSB.hpp"
#include "AMPLModel.hpp"

Level Logger::logger_level = INFO;

// evaluate the gradient of the augmented Lagrangian at x,y, where y=constraint_multipliers
std::vector<double> compute_augmented_lagrangian_gradient(Problem& problem, std::map<int,int>& slacked_constraints, std::vector<double>& x, std::vector<double>& constraints, std::vector<double>& constraint_multipliers, double penalty_parameter) {
    // start with gradient of the objective
    std::vector<double> augmented_lagrangian_gradient = problem.objective_dense_gradient(x);
    // gradient of the constraints wrt the variables
    for (int j = 0; j < problem.number_constraints; j++) {
        double constraint_value;
        try {
            // inequality constraint: need to subtract slacks
            int current_slack = slacked_constraints[j];
            constraint_value = constraints[j] - x[problem.number_variables + current_slack];
        }
        catch (std::out_of_range) {
            // equality constraint
            constraint_value = constraints[j] - problem.constraint_lb[j];
        }
        double factor = constraint_multipliers[j] - penalty_parameter*constraint_value;
        // add the gradient contribution from the constraints
        std::vector<double> constraint_gradient = problem.constraint_dense_gradient(j, x);
        for (int i = 0; i < problem.number_variables; i++) {
            augmented_lagrangian_gradient[i] -= factor*constraint_gradient[i];
        }
    }
    // gradient of the constraints wrt the slacks
    for (std::pair<const int, int> element: slacked_constraints) {
        int j = element.first;                // index of the constraint
        int current_slack = element.second;   // index of the slack in [0, number_of_slacks[
        double constraint_value = constraints[j] - x[problem.number_variables + current_slack];
        double derivative = constraint_multipliers[j] - penalty_parameter*constraint_value;
        augmented_lagrangian_gradient.push_back(derivative); // sticks gradient terms at end of n (number of vars) grad.
    }
    return augmented_lagrangian_gradient;
}

// constraint violation
double compute_eta(Problem& problem, std::map<int,int>& slacked_constraints, std::vector<double>& x, std::vector<double>& constraints) {
    double constraint_violation = 0.;
   
    for (int j = 0; j < problem.number_constraints; j++) {
        double constraint_value;
        try {
            // inequality constraint: need to subtract slacks
            int current_slack = slacked_constraints[j];
            constraint_value = constraints[j] - x[problem.number_variables + current_slack];
        }
        catch (std::out_of_range) {
            // equality constraint
            constraint_value = constraints[j] - problem.constraint_lb[j];
        }
        constraint_violation += std::abs(constraint_value);
    }
    return constraint_violation;
}

// residual of first-order conditions
double compute_omega(Problem& problem, std::map<int,int>& slacked_constraints, std::vector<double>& x, std::vector<double>& constraints, std::vector<double>& constraint_multipliers, double penalty_parameter) {
    // compute the AL gradient
    std::vector<double> augmented_lagrangian_gradient = compute_augmented_lagrangian_gradient(problem, slacked_constraints, x, constraints, constraint_multipliers, penalty_parameter);
    
    double residual = 0.;
    for (unsigned int k = 0; k < x.size(); k++) {
        residual += std::abs(std::min(x[k], augmented_lagrangian_gradient[k]));
    }
    return residual;
}

int main(int argc, char* argv[]) {
    // create the problem
    std::string problem_name = std::string(argv[argc - 1]);
    AMPLModel problem = AMPLModel(problem_name);
    
    // initial primal and dual points
    std::vector<double> x = problem.primal_initial_solution();
    std::vector<double> constraint_multipliers = problem.dual_initial_solution();
    std::vector<double> bound_multipliers(problem.number_variables);
    
    // identify the inequality constraint slacks
    std::vector<double> constraints = problem.evaluate_constraints(x);
    std::map<int,int> slacked_constraints;
    int current_slack = 0;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] != EQUAL_BOUNDS) {
            slacked_constraints[j] = current_slack;
            current_slack++;
            // add slack as a primal variable
            x.push_back(constraints[j]);
        }
    }
    
    // initialize the AL
    double penalty_parameter = 1.;
    // filter
    FilterConstants filter_constants = {0.999, 0.001}; // beta and gamma
    Filter filter(filter_constants);
    // double upper_bound = std::max(this->constants.ubd, this->constants.fact * first_iterate.feasibility_measure);
    // filter.upper_bound = upper_bound;
    // initialize the filter with initial point
    double eta_0 = compute_eta(problem, slacked_constraints, x, constraints);
    double omega_0 = compute_omega(problem, slacked_constraints, x, constraints, constraint_multipliers, penalty_parameter);
    filter.add(eta_0, omega_0);
    std::cout << "Filter entries: " << eta_0 << " " << omega_0 << "\n";
    
    // create the NLP solver
    int limited_memory_size = 5;
    LBFGSB solver(limited_memory_size);
    
    // optimization loop    
    bool optimal = false;
    while (!optimal) {
        bool restoration_phase = false;
        bool filter_acceptable = false;
        
        double eta = 0.;
        double omega = 0.;
        while (!filter_acceptable) {
            // approximately minimize Augmented Lagrangian subproblem
            // TODO
            std::vector<double> trial_x = x;
            
            if (false) { // restoration switching condition (3.14) or (3.15) holds
                restoration_phase = true;
                penalty_parameter *= 2.;
                // Switch to restoration phase to find filter-acceptable point
            }
            else {
                // provisionally update multipliers
                std::vector<double> trial_constraints = problem.evaluate_constraints(trial_x);
                std::vector<double> trial_constraint_multipliers(constraint_multipliers.size());
                for (unsigned int j = 0; j < constraint_multipliers.size(); j++) {
                    trial_constraint_multipliers[j] = constraint_multipliers[j] - penalty_parameter*trial_constraints[j];
                }
                
                // compute filter entries
                eta = compute_eta(problem, slacked_constraints, trial_x, trial_constraints);
                omega = compute_omega(problem, slacked_constraints, trial_x, trial_constraints, trial_constraint_multipliers, penalty_parameter);
                
                // test filter acceptance
                filter_acceptable = true;
            }
        }
        // optional: perform second-order correction + LS + recompute eta and omega
        
        if (eta > 0.) {
            // add entries to filter
            filter.add(eta, omega);
        }
        // optional: penalty update
        double penalty_threshold = 1.;
        if (!restoration_phase && penalty_parameter < penalty_threshold) {
            // increase penalty
        }
        // test optimality
        optimal = true;
    }    
    
    return 0;
}
