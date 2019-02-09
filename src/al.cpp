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
        FilterAugmentedLagrangian();
        
        double compute_eta(std::vector<double>& x, std::vector<double>& constraints);
        double compute_omega(Problem& problem, std::vector<double>& x, std::vector<double>& gradient);
        std::vector<double> compute_constraint_multipliers(Problem& problem, std::vector<double>& constraints, std::vector<double>& constraint_multipliers);
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
        f -= constraint_multipliers[j]*constraints[j]; // f = f - lambda[j]*c[j]
        f += penalty_parameter/2.*constraints[j]*constraints[j]; // f = f + rho/2*(c[j])^2 = augmented Lagrangian
    }
    return f;
}

// evaluate the gradient of the augmented Lagrangian at x,y, where y=constraint_multipliers
std::vector<double> compute_augmented_lagrangian_gradient(Problem& problem, std::vector<double>& x, std::vector<double>& constraints, std::vector<double>& constraint_multipliers, double penalty_parameter) {
    // gradient of the objective
    std::vector<double> augmented_lagrangian_gradient = problem.objective_dense_gradient(x);
    
    // gradient of the constraints wrt original variables
    for (unsigned int j = 0; j < constraints.size(); j++) {
        // dense gradient
        std::vector<double> constraint_gradient = problem.constraint_dense_gradient(j, x);
        double factor = constraint_multipliers[j] - penalty_parameter*constraints[j];
        for (int i = 0; i < problem.number_variables; i++) {
            augmented_lagrangian_gradient[i] -= factor*constraint_gradient[i];
        }
    }
    // gradient of the inequality constraints wrt the slacks
    for (std::pair<const int, int> element: problem.inequality_constraints) {
        int j = element.first; // index of the constraint
        double derivative = constraint_multipliers[j] - penalty_parameter*constraints[j];
        augmented_lagrangian_gradient.push_back(derivative); // sticks gradient terms at end of n (number of vars) grad.
    }
    return augmented_lagrangian_gradient;
}

// constraint violation (infeasibility)
double FilterAugmentedLagrangian::compute_eta(std::vector<double>& x, std::vector<double>& constraints) {
    // compute the constraint violation
    double constraint_violation = 0.;
    for (unsigned int j = 0; j < constraints.size(); j++) {
        constraint_violation += std::abs(constraints[j]);
    }
    return constraint_violation;
}

// residual of first-order conditions
double FilterAugmentedLagrangian::compute_omega(Problem& problem, std::vector<double>& x, std::vector<double>& gradient) {
    // compute the residual
    double residual = 0.;
    for (int i = 0; i < problem.number_variables; i++) {
        residual += std::abs(std::min(x[i] - problem.variable_lb[i], std::max(x[i] - problem.variable_ub[i], gradient[i])));
    }
    for (std::pair<const int, int> element: problem.inequality_constraints) {
        int j = element.first; // index of the constraint
        int slack_index = element.second;
        double slack_value = x[problem.number_variables + slack_index];
        residual += std::abs(std::min(slack_value - problem.constraint_lb[j], std::max(slack_value - problem.constraint_ub[j], gradient[problem.number_variables + slack_index])));
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

std::vector<double> FilterAugmentedLagrangian::compute_constraint_multipliers(Problem& problem, std::vector<double>& constraints, std::vector<double>& constraint_multipliers) {
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
    std::vector<double> constraint_multipliers = problem.dual_initial_solution();
    
    // add slacks as primal variables
    std::vector<double> original_constraints = problem.evaluate_constraints(x);
    for (std::pair<const int, int> element: problem.inequality_constraints) {
        int j = element.first; // index of the constraint
        x.push_back(original_constraints[j]);
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
    
    // compute the AL gradient
    std::vector<double> augmented_lagrangian_gradient = compute_augmented_lagrangian_gradient(problem, x, constraints, constraint_multipliers, this->penalty_parameter); // gradient of f wrt primal variables
    std::cout << "Initial AL gradient: "; print_vector(std::cout, augmented_lagrangian_gradient);
    double eta_0 = this->compute_eta(x, constraints);
    double omega_0 = this->compute_omega(problem, x, augmented_lagrangian_gradient);
    filter.add(eta_0, omega_0);
    std::cout << "Initial filter entries: " << eta_0 << " " << omega_0 << "\n";
    double upper_bound = std::max(100., 1.25*eta_0);
    filter.upper_bound = upper_bound;

    // optimization loop
    double sigma = 0.1; // sufficient reduction
    bool optimal = false;
    int iterations = 0;
    while (!optimal && iterations < 10000) {
        bool restoration_phase = false;

        double eta = 0.;
        double omega = 0.;
        bool filter_acceptable = false;
        strcpy(this->task_, "START");
        while (!filter_acceptable) {
            std::cout << "\n## Starting BFGS from "; print_vector(std::cout, x);
            std::cout << "Filter upper bound is " << filter.upper_bound << "\n";
            std::cout << "Penalty parameter is " << this->penalty_parameter << "\n";
            std::cout << "Current constraints multipliers: "; print_vector(std::cout, constraint_multipliers);
            
            std::vector<double> constraints = compute_constraints(problem, x);
            double initial_augmented_lagrangian = compute_augmented_lagrangian(problem, x, constraints, constraint_multipliers, this->penalty_parameter);
             // compute the AL gradient
            std::vector<double> initial_augmented_lagrangian_gradient = compute_augmented_lagrangian_gradient(problem, x, constraints, constraint_multipliers, this->penalty_parameter);
            double initial_omega = compute_omega(problem, x, initial_augmented_lagrangian_gradient);
            
            /* memory allocation for L-BFGS-B (limited memory & factors allocation) */
            double augmented_lagrangian; // objective
            std::vector<double> wa(this->limited_memory_size*(2*n + 11*this->limited_memory_size + 8) + 5*n);
            std::vector<int> iwa(3*n);
            /************************start BFGS**********************/
            // approximately minimize Augmented Lagrangian subproblem
            //strcpy(this->task_, "START");
            std::vector<double> trial_x(x);
            bool stop_bfgs = false;
            int bfgs_iteration = 0;
            
            // optimization loop (lbfgsb.f uses reverse communication to get function and gradient values)
            while (!stop_bfgs) {
                /* call L-BFGS-B */
                setulb_(&n, &this->limited_memory_size, trial_x.data(), l.data(), u.data(), nbd.data(), &augmented_lagrangian, augmented_lagrangian_gradient.data(), &this->factr_, &this->pgtol_, wa.data(), iwa.data(), this->task_, &this->iprint_, this->csave_, this->lsave_, this->isave_, this->dsave_);
                
                // evaluate Augmented Lagrangian and its gradient
                if (strncmp(this->task_, "FG", 2) == 0) {
                    std::cout << "x: "; print_vector(std::cout, trial_x);
                    std::vector<double> trial_constraints = compute_constraints(problem, trial_x);
                    //std::cout << "constraints: "; print_vector(std::cout, trial_constraints);
                    augmented_lagrangian = compute_augmented_lagrangian(problem, trial_x, trial_constraints, constraint_multipliers, this->penalty_parameter);
                    augmented_lagrangian_gradient = compute_augmented_lagrangian_gradient(problem, trial_x, trial_constraints, constraint_multipliers, this->penalty_parameter);
                    std::cout << "f is " << augmented_lagrangian << "\n";
                    std::cout << "g is "; print_vector(std::cout, augmented_lagrangian_gradient);
                    
                    // termination test
                    if (0 < bfgs_iteration) {
                        // sufficient reduction (3.16)
                        if (initial_augmented_lagrangian - augmented_lagrangian >= sigma*initial_omega) {
                            std::cout << "BFGS termination: sufficient reduction of AL\n";
                            stop_bfgs = true;
                        }
                    }
                    bfgs_iteration++;
                }
                stop_bfgs = stop_bfgs || !(strncmp(this->task_, "FG", 2) == 0 || strncmp(this->task_, "NEW_X", 5) == 0 || strncmp(this->task_, "START", 5) == 0);
            }
            /************************end BFGS**********************/
            std::cout << "BFGS exited with solution x opt: "; print_vector(std::cout, trial_x);
            std::vector<double> trial_constraints = compute_constraints(problem, trial_x);
            std::cout << "Trial constraints: "; print_vector(std::cout, trial_constraints);
            
            // compute eta
            eta = compute_eta(trial_x, trial_constraints);
            
            // test restoration switching conditions
            bool switch_to_restoration = (eta >= beta*filter.upper_bound);
            if (!switch_to_restoration) {
                // second condition
                std::vector<std::map<int, double> > constraints_jacobian = problem.constraints_sparse_jacobian(trial_x);
                std::vector<double> restoration_gradient(trial_x.size());
                for (int j = 0; j < problem.number_constraints; j++) {
                    for (std::pair<const int, double> element: constraints_jacobian[j]) {
                        int i = element.first; // index of the variable
                        double derivative = element.second;
                        restoration_gradient[i] += 2.*derivative*trial_constraints[j];
                    }
                    try {
                        // inequality constraint: need to subtract slack values
                        int slack_index = problem.inequality_constraints.at(j);
                        restoration_gradient[problem.number_variables + slack_index] = -2.*trial_constraints[j];
                    }
                    catch (std::out_of_range) {}
                }
                std::cout << "Restoration gradient: "; print_vector(std::cout, restoration_gradient);
                double omega_restoration = this->compute_omega(problem, trial_x, restoration_gradient);
                std::cout << "Restoration omega = " << omega_restoration << "\n";
                std::cout << "Constraint violation eta = " << eta << "\n";
                switch_to_restoration = (omega_restoration < 1e-8 && eta >= beta*filter.infeasibility_measures[0]);
            }
            
            if (switch_to_restoration) { // restoration switching condition (3.14) or (3.15) holds
                restoration_phase = true;
                this->penalty_parameter *= 2.;
                std::cout << "eta_min was " << filter.infeasibility_measures[0] << "\n";
                // switch to restoration phase to find filter-acceptable point
                throw std::logic_error("Switching to restoration phase");
            }
            else {
                // compute trial multipliers
                std::vector<double> trial_constraint_multipliers = this->compute_constraint_multipliers(problem, trial_constraints, constraint_multipliers);
                std::cout << "Trial constraint multipliers: "; print_vector(std::cout, trial_constraint_multipliers);

                // compute omega
                std::vector<double> trial_augmented_lagrangian_gradient = compute_augmented_lagrangian_gradient(problem, trial_x, trial_constraints, trial_constraint_multipliers, this->penalty_parameter);
                omega = compute_omega(problem, trial_x, trial_augmented_lagrangian_gradient);
                std::cout << "Filter entries: " << eta << " " << omega << "\n";
                
                // test filter acceptance
                filter_acceptable = filter.query(eta, omega);
                std::cout << "Filter acceptable? " << filter_acceptable << "\n";
                // update
                x = trial_x;
                constraint_multipliers = trial_constraint_multipliers;
            }
        }
        // optional: perform second-order correction + LS + recompute eta and omega
        
        if (eta > 0.) {
            // add entries to filter
            std::cout << "Adding entries to the filter\n";
            filter.add(eta, omega);
        }
        // optional: penalty update
        double penalty_threshold = 1.;
        if (!restoration_phase && this->penalty_parameter < penalty_threshold) {
            // increase penalty
            this->penalty_parameter = 2.*penalty_threshold;
        }
        // TODO optimality test
        iterations++;
        if (eta <= 1e-6 && omega <= 1e-5) {
            std::cout << "x is optimal: "; print_vector(std::cout, x);
            optimal = true;
        }
    }
    
    std::cout << "\nSummary:\n";
    std::cout << "Objective evaluations: " << problem.number_eval_objective << "\n";
    std::cout << "Constraint evaluations: " << problem.number_eval_constraints << "\n";
    
    return 0;
}

int main(int argc, char* argv[]) {
    std::string problem_name = std::string(argv[argc - 1]);
    FilterAugmentedLagrangian filter_al;
    return filter_al.solve(problem_name);
}
