#include <iostream>
#include <vector>
#include "Logger.hpp"
#include "Filter.hpp"
#include "AMPLModel.hpp"
#include "LocalSolution.hpp"

Level Logger::logger_level = INFO;

// fortran interface to L-BFGS-B
extern "C" {
    void setulb_(int *n, int *m, double *x, double *l, double *u, int *nbd, double *f, double *g,
            double *factr, double *pgtol, double *wa, int *iwa, char *task, int *iprint,
            char *csave, int *lsave, int *isave, double *dsave);
}

class FilterAugmentedLagrangian {
    public:
        std::map<int,int> inequality_constraints;
        double penalty_parameter;
        FilterAugmentedLagrangian();
        
        double compute_eta(std::vector<double>& x, std::vector<double>& constraints);
        double compute_omega(Problem& problem, std::vector<double>& x, std::vector<double>& constraints, std::vector<double>& constraint_multipliers, std::vector<double>& bound_multipliers);
        std::vector<double> compute_constraint_multipliers(Problem& problem, std::vector<double>& x, std::vector<double>& constraints, std::vector<double>& constraint_multipliers);
        int solve(std::string problem_name);
    
    private:
        int limited_memory_size;
        /* Fortran parameters needed by lbfgsb.f */
        char task_[60];
        char csave_[60];
        int lsave_[4];
        int isave_[44];
        double dsave_[29];
        int iprint_ = -1;
        double factr_ = 1e5;
        double pgtol_ = 1e-5;
};

FilterAugmentedLagrangian::FilterAugmentedLagrangian(): penalty_parameter(200.), limited_memory_size(5) {
}

// evaluate the augmented Lagrangian at x,y, where y=constraint_multipliers
double compute_augmented_lagrangian(Problem& problem, std::vector<double>& x, std::vector<double>& constraints, std::vector<double>& constraint_multipliers, double penalty_parameter) {
    // contribution of the objective
    double f = problem.objective(x);
    
    // contribution of the constraints
    for (unsigned int j = 0; j < constraints.size(); j++) {
        f -= constraint_multipliers[j]*constraints[j]; // f = f - lambda[i]*c[i]
        f += penalty_parameter/2.*constraints[j]*constraints[j]; // f = f + rho/2*(c[i])^2 = augmented Lagrangian
    }
    return f;
}

double compute_augmented_lagrangian2(double (*original_objective)(std::vector<double>&), std::vector<double>& x, std::vector<double>& constraints, std::vector<double>& constraint_multipliers, double penalty_parameter) {
    // contribution of the objective
    double f = original_objective(x);
    
    // contribution of the constraints
    for (unsigned int j = 0; j < constraints.size(); j++) {
        f -= constraint_multipliers[j]*constraints[j]; // f = f - lambda[i]*c[i]
        f += penalty_parameter/2.*constraints[j]*constraints[j]; // f = f + rho/2*(c[i])^2 = augmented Lagrangian
    }
    return f;
}

// evaluate the gradient of the augmented Lagrangian at x,y, where y=constraint_multipliers
std::vector<double> compute_augmented_lagrangian_gradient(Problem& problem, std::map<int,int>& inequality_constraints, std::vector<double>& x, std::vector<double>& constraints, std::vector<double>& constraint_multipliers, double penalty_parameter) {
    // start with gradient of the objective
    std::vector<double> augmented_lagrangian_gradient = problem.objective_dense_gradient(x);
    // gradient of the constraints wrt the variables
    for (unsigned int j = 0; j < constraints.size(); j++) {
        double factor = constraint_multipliers[j] - penalty_parameter*constraints[j];
        // gradient contribution from the constraints
        std::vector<double> constraint_gradient = problem.constraint_dense_gradient(j, x);
        for (int i = 0; i < problem.number_variables; i++) {
            augmented_lagrangian_gradient[i] -= factor*constraint_gradient[i];
        }
    }
    // gradient of the constraints wrt the slacks
    for (std::pair<const int, int> element: inequality_constraints) {
        int j = element.first; // index of the constraint
        double derivative = constraint_multipliers[j] - penalty_parameter*constraints[j];
        augmented_lagrangian_gradient.push_back(derivative); // sticks gradient terms at end of n (number of vars) grad.
    }
    return augmented_lagrangian_gradient;
}

// constraint violation
double FilterAugmentedLagrangian::compute_eta(std::vector<double>& x, std::vector<double>& constraints) {
    double constraint_violation = 0.;
   
    for (unsigned int j = 0; j < constraints.size(); j++) {
        constraint_violation += std::abs(constraints[j]);
    }
    return constraint_violation;
}

// residual of first-order conditions
double FilterAugmentedLagrangian::compute_omega(Problem& problem, std::vector<double>& x, std::vector<double>& constraints, std::vector<double>& constraint_multipliers, std::vector<double>& bound_multipliers) {
    // compute the AL gradient
    std::vector<double> augmented_lagrangian_gradient = compute_augmented_lagrangian_gradient(problem, this->inequality_constraints, x, constraints, constraint_multipliers, this->penalty_parameter);
    
    // multipliers of the bound constraints
    for (int i = 0; i < problem.number_variables; i++) {
        augmented_lagrangian_gradient[i] -= bound_multipliers[i];
    }
    for (std::pair<const int, int> element: inequality_constraints) {
        //int j = element.first; // index of the constraint
        int current_slack = element.second;
        augmented_lagrangian_gradient[problem.number_variables + current_slack] -= bound_multipliers[problem.number_variables + current_slack];
    }
    
    // compute the residual
    double residual = 0.;
    for (int i = 0; i < problem.number_variables; i++) {
        residual += std::abs(std::min(problem.variable_ub[i] - x[i], std::min(x[i] - problem.variable_lb[i], augmented_lagrangian_gradient[i])));
    }
    for (std::pair<const int, int> element: inequality_constraints) {
        int j = element.first; // index of the constraint
        int current_slack = element.second;
        double slack_value = x[problem.number_variables + current_slack];
        residual += std::abs(std::min(problem.constraint_ub[j] - slack_value, std::min(slack_value - problem.constraint_lb[j], augmented_lagrangian_gradient[problem.number_variables + current_slack])));
    }
    return residual;
}

std::vector<double> compute_constraints(Problem& problem, std::map<int,int>& inequality_constraints, std::vector<double>& x) {
    std::vector<double> original_constraints = problem.evaluate_constraints(x);
    std::vector<double> constraints(problem.number_constraints);
    for (int j = 0; j < problem.number_constraints; j++) {
        try {
            // inequality constraint: need to subtract slack values
            int current_slack = inequality_constraints[j];
            constraints[j] = original_constraints[j] - x[problem.number_variables + current_slack];
        }
        catch (std::out_of_range) {
            // equality constraint
            constraints[j] = original_constraints[j] - problem.constraint_lb[j];
        }
    }
    return constraints;
}

std::vector<double> FilterAugmentedLagrangian::compute_constraint_multipliers(Problem& problem, std::vector<double>& x, std::vector<double>& constraints, std::vector<double>& constraint_multipliers) {
    /* compute the reduced gradient */
    //double reduced_gradient = 0.;
    //for (int i = 0; i < x.size(); i++) {
    //    std::cout << x[i] << " in [" << l[i] << ", " << u[i] << "]\tmultiplier: " << g[i] << "\n";
    //    reduced_gradient += std::abs(std::min(x[i] - l[i], u[i] - x[i]) * g[i]);
    //}; // end for
    //std::cout << "Reduced Gradient Norm = " << reduced_gradient << "\n";
    
    /* compute the new multipliers by using first-order update form: y_trial = y - rho*c */
    std::vector<double> trial_constraint_multipliers(constraint_multipliers);
    for (int j = 0; j < problem.number_constraints; j++) {
        trial_constraint_multipliers[j] -= this->penalty_parameter*constraints[j];
    }
    return trial_constraint_multipliers;
}

int FilterAugmentedLagrangian::solve(std::string problem_name) {
    // create the problem
    AMPLModel problem = AMPLModel(problem_name);
    
    // initial primal and dual points
    std::vector<double> x = problem.primal_initial_solution();
    std::vector<double> bound_multipliers(problem.number_variables);
    std::vector<double> constraint_multipliers = problem.dual_initial_solution();
    
    // identify the inequality constraint slacks
    std::vector<double> original_constraints = problem.evaluate_constraints(x);
    int current_slack = 0;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] != EQUAL_BOUNDS) {
            this->inequality_constraints[j] = current_slack;
            current_slack++;
            // add slack as a primal variable
            x.push_back(original_constraints[j]);
            bound_multipliers.push_back(0.);
        }
    }

    // compute bounds and variable status
    int n = x.size();
    std::vector<double> l(n);
    std::vector<double> u(n);
    std::vector<ConstraintType> variable_status(n);
    for (int i = 0; i < problem.number_variables; i++) {
        l[i] = problem.variable_lb[i];
        u[i] = problem.variable_ub[i];
        variable_status[i] = problem.variable_status[i];
    }
    for (std::pair<const int, int> element: this->inequality_constraints) {
        int j = element.first;
        int current_slack = element.second;
        l[problem.number_variables + current_slack] = problem.constraint_lb[j];
        u[problem.number_variables + current_slack] = problem.constraint_ub[j];
        variable_status[problem.number_variables + current_slack] = problem.constraint_status[j];
    }
    // set the type of bounds
    std::vector<int> nbd(n);
    for (int i = 0; i < n; i++) {
        if (variable_status[i] == UNBOUNDED) {
            nbd[i] = 0;
        }
        else if (variable_status[i] == BOUNDED_LOWER) {
            nbd[i] = 1;
        }
        else if (variable_status[i] == BOUNDED_UPPER) {
            nbd[i] = 3;
        }
        else { // two bounds
            nbd[i] = 2;
        }
    }
    
    // initialize the AL filter
    FilterConstants filter_constants = {0.999, 0.001}; // beta and gamma
    Filter filter(filter_constants);
    // double upper_bound = std::max(this->constants.ubd, this->constants.fact * first_iterate.feasibility_measure);
    // filter.upper_bound = upper_bound;
    // initialize the filter with initial point's entries
    std::vector<double> constraints = compute_constraints(problem, this->inequality_constraints, x);
    double eta_0 = this->compute_eta(x, constraints);
    double omega_0 = this->compute_omega(problem, x, constraints, constraint_multipliers, bound_multipliers);
    filter.add(eta_0, omega_0);
    std::cout << "Initial filter entries: " << eta_0 << " " << omega_0 << "\n";
    
    // create the NLP solver
    /* memory allocation for L-BFGS-B (limited memory & factors allocation) */
    std::vector<double> wa(this->limited_memory_size*(2*n + 11*this->limited_memory_size + 8) + 5*n);
    std::vector<int> iwa(3*n);
    double f; // objective
    std::vector<double> g(n); // gradient of f wrt primal variables

    // optimization loop    
    bool optimal = false;
    int iterations = 0;
    while (!optimal) {
        bool restoration_phase = false;
        bool filter_acceptable = false;
        
        double eta = 0.;
        double omega = 0.;
        while (!filter_acceptable) {
            // approximately minimize Augmented Lagrangian subproblem
            strcpy(this->task_, "START");
            std::vector<double> x_bfgs(x);
            
            bool stop = false;
            // optimization loop (lbfgsb.f uses reverse communication to get function and gradient values)
            while (!stop) {
                /* call L-BFGS-B */
                setulb_(&n, &this->limited_memory_size, x_bfgs.data(), l.data(), u.data(), nbd.data(), &f, g.data(), &this->factr_, &this->pgtol_, wa.data(), iwa.data(), this->task_, &this->iprint_, this->csave_, this->lsave_, this->isave_, this->dsave_);
                // evaluate Augmented Lagrangian and its gradient
                if (strncmp(this->task_, "FG", 2) == 0) {
                    std::cout << "x: "; print_vector(std::cout, x_bfgs);
                    std::vector<double> constraints_bfgs = compute_constraints(problem, this->inequality_constraints, x_bfgs);
                    f = compute_augmented_lagrangian(problem, x_bfgs, constraints_bfgs, constraint_multipliers, this->penalty_parameter);
                    g = compute_augmented_lagrangian_gradient(problem, this->inequality_constraints, x_bfgs, constraints_bfgs, constraint_multipliers, this->penalty_parameter);
                    std::cout << "f is " << f << "\n";
                    std::cout << "g is "; print_vector(std::cout, g);
                }
                stop = !(strncmp(this->task_, "FG", 2) == 0 || strncmp(this->task_, "NEW_X", 5) == 0 || strncmp(this->task_, "START", 5) == 0);
            }
            
            // no constraint: empty constraint multipliers
            std::vector<double> solution_constraint_multipliers;
            // create local solution from primal and dual variables
            LocalSolution solution(x_bfgs, g, solution_constraint_multipliers);
            solution.status = OPTIMAL;
            
            std::cout << "L-BFGS-B exited with solution\n";
            std::vector<double> trial_x = solution.x;
            std::cout << "Bound multipliers: "; print_vector(std::cout, solution.bound_multipliers);
            
            if (false) { // restoration switching condition (3.14) or (3.15) holds
                restoration_phase = true;
                this->penalty_parameter *= 2.;
                // Switch to restoration phase to find filter-acceptable point
            }
            else {
                std::vector<double> trial_constraints = compute_constraints(problem, this->inequality_constraints, trial_x);
                
                // provisionally update multipliers
                std::vector<double> trial_constraint_multipliers = this->compute_constraint_multipliers(problem, trial_x, trial_constraints, constraint_multipliers);
                std::cout << "Trial constraint multipliers: "; print_vector(std::cout, trial_constraint_multipliers);

                // compute filter entries
                eta = compute_eta(trial_x, trial_constraints);
                omega = compute_omega(problem, trial_x, trial_constraints, trial_constraint_multipliers, solution.bound_multipliers);
                std::cout << "Filter entries: " << eta << " " << omega << "\n";
                
                // test filter acceptance
                filter_acceptable = filter.query(eta, omega);
                std::cout << "Filter acceptable? " << filter_acceptable << "\n";
                // update
                x = trial_x;
                bound_multipliers = solution.bound_multipliers;
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
        // TODO test optimality
        iterations++;
        optimal = true;
        if (eta <= 1e-6 && omega <= 1e-6) {
            std::cout << "x is optimal\n";
            optimal = true;
        }
    }
    
    return 0;
}

int main(int argc, char* argv[]) {
    std::string problem_name = std::string(argv[argc - 1]);
    FilterAugmentedLagrangian filter_al;
    return filter_al.solve(problem_name);
}
