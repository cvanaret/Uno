#include <iostream>
#include <vector>
#include "Logger.hpp"
#include "Filter.hpp"
#include "AMPLModel.hpp"
#include "Utils.hpp"

Level Logger::logger_level = INFO;

// fortran interface to L-BFGS-B
extern "C" {
    void setulb_(int *n, int *m, double *x, double *l, double *u, int *nbd, double *f, double *g,
            double *factr, double *pgtol, double *wa, int *iwa, char *task, int *iprint,
            char *csave, int *lsave, int *isave, double *dsave);
}

class FilterAugmentedLagrangian {
    public:
        double penalty_parameter;
        double penalty_update_factor;
        FilterAugmentedLagrangian();
        
        double compute_eta(std::vector<double>& x, std::vector<double>& constraints);
        double compute_omega(Problem& problem, std::vector<double>& x, std::vector<double>& gradient);
        std::vector<double> compute_constraint_multipliers(Problem& problem, std::vector<double>& trial_constraints, std::vector<double>& constraint_multipliers);
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

FilterAugmentedLagrangian::FilterAugmentedLagrangian(): penalty_parameter(10.), penalty_update_factor(10.), limited_memory_size(5) {
}

// evaluate the augmented Lagrangian at x,y, where y=constraint_multipliers
double compute_augmented_lagrangian(Problem& problem, std::vector<double>& x, std::vector<double>& constraints, std::vector<double>& constraint_multipliers, double penalty_parameter) {
    // contribution of the objective
    double augmented_lagrangian = problem.objective(x);
    
    // contribution of the constraints
    for (unsigned int j = 0; j < constraints.size(); j++) {
        augmented_lagrangian -= constraint_multipliers[j]*constraints[j]; // f = f - lambda[j]*c[j]
        augmented_lagrangian += penalty_parameter/2.*constraints[j]*constraints[j]; // f = f + rho/2 c[j]^2
    }
    return augmented_lagrangian;
}

// evaluate the gradient of the augmented Lagrangian at x,y, where y=constraint_multipliers
std::vector<double> compute_augmented_lagrangian_gradient(Problem& problem, std::vector<double>& x, std::vector<double>& constraints, std::vector<double>& constraint_multipliers, double penalty_parameter) {
    // gradient of the objective
    std::vector<double> augmented_lagrangian_gradient = problem.objective_dense_gradient(x);
    
    // gradient of the constraints wrt original variables
    for (unsigned int j = 0; j < constraints.size(); j++) {
        double factor_j = constraint_multipliers[j] - penalty_parameter*constraints[j];
        // dense gradient (TODO: sparse)
        std::vector<double> constraint_j_gradient = problem.constraint_dense_gradient(j, x);
        for (int i = 0; i < problem.number_variables; i++) {
            augmented_lagrangian_gradient[i] -= factor_j*constraint_j_gradient[i];
        }
    }
    // gradient of the inequality constraints wrt slacks
    for (std::pair<const int, int> element: problem.inequality_constraints) {
        int j = element.first; // index of the constraint
        double derivative = constraint_multipliers[j] - penalty_parameter*constraints[j];
        augmented_lagrangian_gradient.push_back(derivative);
    }
    return augmented_lagrangian_gradient;
}

std::vector<double> compute_lagrangian_gradient(Problem& problem, std::vector<double>& x, std::vector<double>& constraint_multipliers) {
    // gradient of the objective
    std::vector<double> lagrangian_gradient = problem.objective_dense_gradient(x);
    
    // gradient of the constraints wrt original variables
    for (int j = 0; j < problem.number_constraints; j++) {
        // dense gradient (TODO: sparse)
        std::vector<double> constraint_j_gradient = problem.constraint_dense_gradient(j, x);
        for (int i = 0; i < problem.number_variables; i++) {
            lagrangian_gradient[i] -= constraint_multipliers[j]*constraint_j_gradient[i];
        }
    }
    // gradient of the inequality constraints wrt slacks
    for (std::pair<const int, int> element: problem.inequality_constraints) {
        int j = element.first; // index of the constraint
        lagrangian_gradient.push_back(constraint_multipliers[j]);
    }
    return lagrangian_gradient;
}

// constraint violation (infeasibility)
double FilterAugmentedLagrangian::compute_eta(std::vector<double>& x, std::vector<double>& constraints) {
    // compute the constraint violation
    //return norm_1(constraints);
    return norm_2_squared(constraints);
}

// residual of first-order conditions
double FilterAugmentedLagrangian::compute_omega(Problem& problem, std::vector<double>& x, std::vector<double>& lagrangian_gradient) {
    // compute the residual
    double residual = 0.;
    for (int i = 0; i < problem.number_variables; i++) {
        double error = std::min(x[i] - problem.variable_lb[i], std::max(x[i] - problem.variable_ub[i], lagrangian_gradient[i]));
        residual += error*error;
    }
    for (std::pair<const int, int> element: problem.inequality_constraints) {
        int j = element.first; // index of the constraint
        int slack_index = element.second;
        double slack_value = x[problem.number_variables + slack_index];
        double error = std::min(slack_value - problem.constraint_lb[j], std::max(slack_value - problem.constraint_ub[j], lagrangian_gradient[problem.number_variables + slack_index]));
        residual += error*error;
    }
    return residual;
}

std::vector<double> compute_constraints(Problem& problem, std::vector<double>& x) {
    std::vector<double> original_constraints = problem.evaluate_constraints(x);
    std::vector<double> constraints(problem.number_constraints);
    for (int j = 0; j < problem.number_constraints; j++) {
        try {
            // inequality constraint: need to subtract slack value
            int slack_index = problem.inequality_constraints.at(j);
            constraints[j] = original_constraints[j] - x[problem.number_variables + slack_index];
        }
        catch (std::out_of_range) {
            // equality constraint
            constraints[j] = original_constraints[j] - problem.constraint_lb[j];
        }
    }
    return constraints;
}

std::vector<double> FilterAugmentedLagrangian::compute_constraint_multipliers(Problem& problem, std::vector<double>& trial_constraints, std::vector<double>& constraint_multipliers) {
    /* compute the new multipliers by using first-order update form: y_trial = y - rho*c */
    std::vector<double> trial_constraint_multipliers(constraint_multipliers);
    for (int j = 0; j < problem.number_constraints; j++) {
        trial_constraint_multipliers[j] -= this->penalty_parameter*trial_constraints[j];
    }
    return trial_constraint_multipliers;
}

int FilterAugmentedLagrangian::solve(std::string problem_name) {
    // create the problem
    AMPLModel problem = AMPLModel(problem_name);
    double epsilon = 1e-7;
    int number_outer_iterations = 1000;

    // initial primal and dual points
    std::vector<double> x = problem.primal_initial_solution();
    std::vector<double> constraint_multipliers = problem.dual_initial_solution();
    
    // add slacks as primal variables
    for (std::pair<const int, int> element: problem.inequality_constraints) {
        int j = element.first; // index of the constraint
        x.push_back(problem.evaluate_constraint(j, x));
    }
    std::vector<double> bound_multipliers = std::vector<double>(x.size());

    // compute bounds and variable status
    int n = x.size();
    std::cout << "n is " << n << "\n";
    std::vector<double> l(n);
    std::vector<double> u(n);
    std::vector<ConstraintType> variable_status(n);
    for (int i = 0; i < problem.number_variables; i++) {
        l[i] = problem.variable_lb[i];
        u[i] = problem.variable_ub[i];
        variable_status[i] = problem.variable_status[i];
    }
    for (std::pair<const int, int> element: problem.inequality_constraints) {
        int j = element.first;
        int slack_index = element.second;
        l[problem.number_variables + slack_index] = problem.constraint_lb[j];
        u[problem.number_variables + slack_index] = problem.constraint_ub[j];
        variable_status[problem.number_variables + slack_index] = problem.constraint_status[j];
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
    double beta = 0.999;
    double gamma = 0.001;
    FilterConstants filter_constants = {beta, gamma};
    Filter filter(filter_constants);
    
    // initialize the filter with initial point's entries
    std::cout << "Initial x: "; print_vector(std::cout, x);
    std::vector<double> constraints = compute_constraints(problem, x);
    std::cout << "Initial constraints: "; print_vector(std::cout, constraints);
    
    // compute the gradient of the augmented Lagrangian wrt primal variables
    std::vector<double> lagrangian_gradient = compute_lagrangian_gradient(problem, x, constraint_multipliers);
    double omega = this->compute_omega(problem, x, lagrangian_gradient);
    double eta = this->compute_eta(x, constraints);
    if (0. < eta) {
        filter.add(eta, omega);
        std::cout << "Initial filter entries: " << eta << " " << omega << "\n";
    }
    filter.upper_bound = std::max(100., 1.25*eta);
    
    /* memory allocation for L-BFGS-B (limited memory & factors allocation) */
    std::vector<double> wa(this->limited_memory_size*(2*n + 11*this->limited_memory_size + 8) + 5*n);
    std::vector<int> iwa(3*n);

    // optimization loop
    double sigma = 0.1; // sufficient reduction
    bool termination = false;
    int iterations = 0;
    strcpy(this->task_, "START");
    
    while (!termination && iterations < number_outer_iterations) {
        bool filter_acceptable = false;
        while (!filter_acceptable && !termination) {
            std::cout << "\n## Starting BFGS from "; print_vector(std::cout, x);
            std::cout << "Filter upper bound is " << filter.upper_bound << "\n";
            std::cout << "Penalty parameter is " << this->penalty_parameter << "\n";
            std::cout << "Current constraints multipliers: "; print_vector(std::cout, constraint_multipliers);
            
            std::vector<double> constraints = compute_constraints(problem, x);
            double augmented_lagrangian = compute_augmented_lagrangian(problem, x, constraints, constraint_multipliers, this->penalty_parameter);
            std::vector<double> augmented_lagrangian_gradient(x.size());
            
            /************************start BFGS**********************/
            // approximately minimize Augmented Lagrangian subproblem
            std::vector<double> trial_x(x);
            bool stop_bfgs = false;
            int bfgs_iteration = 0;
            double initial_augmented_lagrangian = augmented_lagrangian;
            
            // BFGS loop (lbfgsb.f uses reverse communication to get function and gradient values)
            while (!stop_bfgs) {
                /* call L-BFGS-B */
                setulb_(&n, &this->limited_memory_size, trial_x.data(), l.data(), u.data(), nbd.data(), &augmented_lagrangian, augmented_lagrangian_gradient.data(), &this->factr_, &this->pgtol_, wa.data(), iwa.data(), this->task_, &this->iprint_, this->csave_, this->lsave_, this->isave_, this->dsave_);
                
                // evaluate Augmented Lagrangian and its gradient
                if (strncmp(this->task_, "FG", 2) == 0) {
                    std::cout << "x: "; print_vector(std::cout, trial_x);
                    std::vector<double> trial_constraints = compute_constraints(problem, trial_x);
                    // compute augmented Lagrangian and its gradient at the point (x_trial, y)
                    if (0 < bfgs_iteration) {
                        augmented_lagrangian = compute_augmented_lagrangian(problem, trial_x, trial_constraints, constraint_multipliers, this->penalty_parameter);
                    }
                    augmented_lagrangian_gradient = compute_augmented_lagrangian_gradient(problem, trial_x, trial_constraints, constraint_multipliers, this->penalty_parameter);
                    std::cout << "f is " << augmented_lagrangian << "\n";
                    std::cout << "g is "; print_vector(std::cout, augmented_lagrangian_gradient);
                    
                    // termination test
                    if (0 < bfgs_iteration) {
                        // sufficient reduction (3.16)
                        if (initial_augmented_lagrangian - augmented_lagrangian >= sigma*omega) {
                            std::cout << "BFGS termination: sufficient reduction of AL between " << initial_augmented_lagrangian << " and " << augmented_lagrangian << "\n";
                            stop_bfgs = true;
                        }
                    }
                    bfgs_iteration++;
                }
                if (strncmp(this->task_, "ABNO", 4) == 0 || strncmp(this->task_, "ERROR", 5) == 0) {
                    std::cout << "Task: " << this->task_ << "\n";
                    throw std::logic_error("Error in BFGS");
                }
                if (strncmp(this->task_, "CONV", 4) == 0) {
                    std::cout << "BFGS converged with stationary point\n";
                    stop_bfgs = true;
                }
            }
            std::cout << "BFGS exited with solution x*: "; print_vector(std::cout, trial_x);
            /************************end BFGS**********************/
            std::vector<double> trial_constraints = compute_constraints(problem, trial_x);
            std::cout << "c(x*) = "; print_vector(std::cout, trial_constraints);
            // compute gradient of bound-constrained augmented Lagrangian at (x*, y)
            augmented_lagrangian_gradient = compute_augmented_lagrangian_gradient(problem, trial_x, trial_constraints, constraint_multipliers, this->penalty_parameter);
            // compute filter entries
            eta = compute_eta(trial_x, trial_constraints);
            omega = compute_omega(problem, trial_x, augmented_lagrangian_gradient);
            std::cout << "Filter entries: (η = " << eta << ", ω = " << omega << ")\n";
            
            // termination test
            if (eta <= epsilon && (strncmp(this->task_, "CONV", 4) == 0 || omega <= epsilon)) {
                // convergence towards a feasible point: check if dual infeasibility is <= epsilon
                if (epsilon < omega) {
                    std::cout << "Subsolver converged, but couldn't drive dual infeasibility below tolerance\n";
                }
                else {
                    std::cout << "Subsolver converged within given tolerance\n";
                }
                termination = true;
            }
            // restoration switching conditions
            else if ((eta >= beta*filter.upper_bound) || (omega <= epsilon && filter.size > 0 && eta >= beta*std::max(filter.infeasibility_measures[0], epsilon))) { 
                // restoration switching condition (3.14) or (3.15) holds
                // switch to restoration phase to find filter-acceptable point
                //throw std::logic_error("Switching to restoration phase");
                std::cout << "  *** Switching to restoration phase\n";
                
                /************************start BFGS**********************/
                // approximately minimize Augmented Lagrangian subproblem
                std::vector<double> trial_x_restoration(trial_x);
                double restoration_objective;
                std::vector<double> restoration_objective_gradient;
                
                bool stop_bfgs_restoration = false;
                strcpy(this->task_, "START");
                
                // BFGS loop (lbfgsb.f uses reverse communication to get function and gradient values)
                while (!stop_bfgs_restoration) {
                    /* call L-BFGS-B */
                    setulb_(&n, &this->limited_memory_size, trial_x_restoration.data(), l.data(), u.data(), nbd.data(), &restoration_objective, restoration_objective_gradient.data(), &this->factr_, &this->pgtol_, wa.data(), iwa.data(), this->task_, &this->iprint_, this->csave_, this->lsave_, this->isave_, this->dsave_);
                    
                    // evaluate Augmented Lagrangian and its gradient
                    if (strncmp(this->task_, "FG", 2) == 0) {
                        std::cout << "x: "; print_vector(std::cout, trial_x_restoration);
                        std::vector<double> constraints = compute_constraints(problem, trial_x_restoration);
                        // compute constraint violation at the point trial_x_restoration
                        restoration_objective = 0.;
                        for (unsigned int j = 0; j < constraints.size(); j++) {
                            restoration_objective += constraints[j]*constraints[j];
                        }
                        // compute constraint violation gradient at the point trial_x_restoration
                        restoration_objective_gradient = std::vector<double>(problem.number_variables);
                        for (unsigned int j = 0; j < constraints.size(); j++) {
                            double factor_j = 2.*constraints[j];
                            std::vector<double> constraint_j_gradient = problem.constraint_dense_gradient(j, trial_x_restoration);
                            for (int i = 0; i < problem.number_variables; i++) {
                                restoration_objective_gradient[i] += factor_j*constraint_j_gradient[i];
                            }
                        }
                        // gradient of the inequality constraints wrt slacks
                        for (std::pair<const int, int> element: problem.inequality_constraints) {
                            int j = element.first; // index of the constraint
                            restoration_objective_gradient.push_back(-2.*constraints[j]);
                        }
                        std::cout << "f is " << restoration_objective << "\n";
                        std::cout << "g is "; print_vector(std::cout, restoration_objective_gradient);
                    }
                    if (strncmp(this->task_, "ABNO", 4) == 0 || strncmp(this->task_, "ERROR", 5) == 0) {
                        std::cout << "Task: " << this->task_ << "\n";
                        throw std::logic_error("Error in BFGS");
                    }
                    if (strncmp(this->task_, "CONV", 4) == 0) {
                        std::cout << "BFGS converged with stationary point\n";
                        stop_bfgs_restoration = true;
                    }
                }
                std::cout << "BFGS exited with solution x*: "; print_vector(std::cout, trial_x_restoration);
                /************************end BFGS**********************/
                
                if (1e-5 < restoration_objective) {
                    std::cout << "Restoration phase ended with non-feasible FJ point\n";
                    termination = true;
                }
                else {
                    std::vector<double> trial_constraints_restoration = compute_constraints(problem, trial_x_restoration);
                    std::cout << "c(x*) = "; print_vector(std::cout, trial_constraints_restoration);
                    // compute gradient of bound-constrained augmented Lagrangian at (x*, y)
                    augmented_lagrangian_gradient = compute_augmented_lagrangian_gradient(problem, trial_x_restoration, trial_constraints_restoration, constraint_multipliers, this->penalty_parameter);
                    // compute filter entries
                    eta = compute_eta(trial_x_restoration, trial_constraints_restoration);
                    omega = compute_omega(problem, trial_x_restoration, augmented_lagrangian_gradient);
                    std::cout << "Restoration filter entries: (η = " << eta << ", ω = " << omega << ")\n";

                    // test filter acceptance
                    filter_acceptable = filter.query(eta, omega);
                    std::cout << "Filter acceptable? " << (filter_acceptable ? "yes" : "no") << "\n";
                    // update
                    x = trial_x_restoration;
                    // TODO check
                    std::vector<double> trial_constraint_multipliers = this->compute_constraint_multipliers(problem, trial_constraints, constraint_multipliers);
                    constraint_multipliers = trial_constraint_multipliers;
                    bound_multipliers = augmented_lagrangian_gradient;
                    this->penalty_parameter *= this->penalty_update_factor;
                    std::cout << "Restoration triggered: increasing parameter\n";
                }
            }
            else {
                // test filter acceptance
                filter_acceptable = filter.query(eta, omega);
                std::cout << "Filter acceptable? " << (filter_acceptable ? "yes" : "no") << "\n";
                // update
                x = trial_x;
                if (filter_acceptable) {
                    // update multipliers
                    std::vector<double> trial_constraint_multipliers = this->compute_constraint_multipliers(problem, trial_constraints, constraint_multipliers);
                    std::cout << "λ(x*) = "; print_vector(std::cout, trial_constraint_multipliers);
                    constraint_multipliers = trial_constraint_multipliers;
                    bound_multipliers = augmented_lagrangian_gradient;
                    strcpy(this->task_, "START");
                }
            }
        }
        
        
        // optional: perform second-order correction + LS + recompute eta and omega
        
        if (0. < eta) {
            std::cout << "Adding entries to the filter\n";
            filter.add(eta, omega);
        }
        iterations++;
    }
    if (iterations >= number_outer_iterations) {
        std::cout << "number_outer_iterations was reached\n";
    }
    
    //if (eta <= epsilon && omega <= epsilon) {
            //std::cout << "\nOPTIMAL SOLUTION\n";
    std::cout << "\nOptimization summary:\n";
    std::cout << "x = "; print_vector(std::cout, x);
    std::cout << "constraint multipliers = "; print_vector(std::cout, constraint_multipliers);
    std::cout << "bound multipliers = "; print_vector(std::cout, bound_multipliers);
    std::cout << "Objective evaluations: " << problem.number_eval_objective << "\n";
    std::cout << "Constraint evaluations: " << problem.number_eval_constraints << "\n";
    
    return 0;
}

int main(int argc, char* argv[]) {
    std::string problem_name = std::string(argv[argc - 1]);
    FilterAugmentedLagrangian filter_al;
    return filter_al.solve(problem_name);
}

// second condition
//std::vector<std::map<int, double> > constraints_jacobian = problem.constraints_sparse_jacobian(trial_x);
//std::vector<double> restoration_gradient(trial_x.size());
//for (int j = 0; j < problem.number_constraints; j++) {
    //double constraint_value = trial_constraints[j];
    //for (std::pair<const int, double> element: constraints_jacobian[j]) {
        //int i = element.first; // index of the variable
        //double derivative = element.second;
        //restoration_gradient[i] += 2.*derivative*constraint_value;
    //}
    //try {
        //// inequality constraint: need to subtract slack values
        //int slack_index = problem.inequality_constraints.at(j);
        //restoration_gradient[problem.number_variables + slack_index] = -2.*constraint_value;
    //}
    //catch (std::out_of_range) {}
//}
//std::cout << "Restoration gradient: "; print_vector(std::cout, restoration_gradient);
//double omega_restoration = this->compute_omega(problem, trial_x, restoration_gradient);
//std::cout << "Restoration omega = " << omega_restoration << "\n";
