#include <iostream>
#include <vector>
#include "Logger.hpp"
#include "Filter.hpp"
#include "LBFGSB.hpp"
#include "AMPLModel.hpp"

Level Logger::logger_level = INFO;

class FilterAugmentedLagrangian {
    public:
        std::map<int,int> slacked_constraints;
        double penalty_parameter;
        FilterAugmentedLagrangian();
        
        // TODO constructor

        double compute_eta(Problem& problem, std::vector<double>& x, std::vector<double>& constraints);
        double compute_omega(Problem& problem, std::vector<double>& x, std::vector<double>& constraints, std::vector<double>& constraint_multipliers);
        int solve(std::string problem_name);
};

FilterAugmentedLagrangian::FilterAugmentedLagrangian(): penalty_parameter(1.) {
}

// evaluate the augmented Lagrangian at x,y, where y=constraint_multipliers
double compute_augmented_lagrangian(Problem& problem, std::map<int,int>& slacked_constraints, std::vector<double>& x, std::vector<double>& constraints, std::vector<double>& constraint_multipliers, double penalty_parameter) {
    // contribution of the objective
    double f = problem.objective(x);
    // contribution of the constraints
    for (int j = 0; j < problem.number_constraints; j++) {
        double constraint_value;
        try {
            // inequality constraint: need to subtract slack values
            int current_slack = slacked_constraints[j];
            constraint_value = constraints[j] - x[problem.number_variables + current_slack];
        }
        catch (std::out_of_range) {
            // equality constraint
            constraint_value = constraints[j] - problem.constraint_lb[j];
        }
        f -= constraint_multipliers[j]*constraint_value;       // f = f - lambda[i]*c[i]
        f += penalty_parameter/2.*constraint_value*constraint_value;   // f = f + rho/2*(c[i])^2 = augmented Lagrangian
    }
    return f;
}

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
double FilterAugmentedLagrangian::compute_eta(Problem& problem, std::vector<double>& x, std::vector<double>& constraints) {
    double constraint_violation = 0.;
   
    for (int j = 0; j < problem.number_constraints; j++) {
        double constraint_value;
        try {
            // inequality constraint: need to subtract slacks
            int current_slack = this->slacked_constraints[j];
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
double FilterAugmentedLagrangian::compute_omega(Problem& problem, std::vector<double>& x, std::vector<double>& constraints, std::vector<double>& constraint_multipliers) {
    // compute the AL gradient
    std::vector<double> augmented_lagrangian_gradient = compute_augmented_lagrangian_gradient(problem, this->slacked_constraints, x, constraints, constraint_multipliers, this->penalty_parameter);
    
    double residual = 0.;
    for (unsigned int k = 0; k < x.size(); k++) {
        residual += std::abs(std::min(x[k], augmented_lagrangian_gradient[k]));
    }
    return residual;
}

int FilterAugmentedLagrangian::solve(std::string problem_name) {
    // create the problem
    AMPLModel problem = AMPLModel(problem_name);
    
    // initial primal and dual points
    std::vector<double> x = problem.primal_initial_solution();
    std::vector<double> bound_multipliers(problem.number_variables);
    std::vector<double> constraint_multipliers = problem.dual_initial_solution();
    Iterate current_iterate(problem, x, bound_multipliers, constraint_multipliers);
    
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
    
    // initialize the AL filter
    FilterConstants filter_constants = {0.999, 0.001}; // beta and gamma
    Filter filter(filter_constants);
    // double upper_bound = std::max(this->constants.ubd, this->constants.fact * first_iterate.feasibility_measure);
    // filter.upper_bound = upper_bound;
    // initialize the filter with initial point's entries
    double eta_0 = this->compute_eta(problem, x, constraints);
    double omega_0 = this->compute_omega(problem, x, constraints, constraint_multipliers);
    filter.add(eta_0, omega_0);
    std::cout << "Filter entries: " << eta_0 << " " << omega_0 << "\n";
    
    // create the NLP solver
    int limited_memory_size = 5;
    LBFGSB nlp_solver(limited_memory_size);
    nlp_solver.initialize(slacked_constraints);
    
    // optimization loop    
    bool optimal = false;
    while (!optimal) {
        bool restoration_phase = false;
        bool filter_acceptable = false;
        
        double eta = 0.;
        double omega = 0.;
        while (!filter_acceptable) {
            // approximately minimize Augmented Lagrangian subproblem
            LocalSolution solution = nlp_solver.solve(problem, current_iterate, compute_augmented_lagrangian, compute_augmented_lagrangian_gradient);
            std::cout << "L-BFGS-B exited with solution\n";
            std::vector<double> trial_x = solution.x;
            
            if (false) { // restoration switching condition (3.14) or (3.15) holds
                restoration_phase = true;
                this->penalty_parameter *= 2.;
                // Switch to restoration phase to find filter-acceptable point
            }
            else {
                // provisionally update multipliers
                std::vector<double> trial_constraints = problem.evaluate_constraints(trial_x);
                std::vector<double> trial_constraint_multipliers(constraint_multipliers.size());
                for (unsigned int j = 0; j < constraint_multipliers.size(); j++) {
                    trial_constraint_multipliers[j] = constraint_multipliers[j] - this->penalty_parameter*trial_constraints[j];
                }
                
                // compute filter entries
                eta = compute_eta(problem, trial_x, trial_constraints);
                omega = compute_omega(problem, trial_x, trial_constraints, trial_constraint_multipliers);
                
                // test filter acceptance
                filter_acceptable = true;
                x = trial_x;
                constraint_multipliers = trial_constraint_multipliers;
            }
        }
        // optional: perform second-order correction + LS + recompute eta and omega
        
        if (eta > 0.) {
            // add entries to filter
            filter.add(eta, omega);
        }
        // optional: penalty update
        double penalty_threshold = 1.;
        if (!restoration_phase && this->penalty_parameter < penalty_threshold) {
            // increase penalty
        }
        // test optimality
        optimal = true;
    }    
    
    return 0;
}

int main(int argc, char* argv[]) {
    std::string problem_name = std::string(argv[argc - 1]);
    FilterAugmentedLagrangian filter_al;
    return filter_al.solve(problem_name);
}
