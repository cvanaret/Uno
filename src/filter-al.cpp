#include <iostream>
#include <vector>
#include <set>
#include "Logger.hpp"
#include "Constraint.hpp"
#include "Filter.hpp"
#include "AMPLModel.hpp"
#include "Utils.hpp"
#include "BQPDSolver.hpp"

Level Logger::logger_level = INFO;

// fortran interface to L-BFGS-B
extern "C" {
    void setulb_(int *n, int *m, double *x, double *l, double *u, int *nbd, double *f, double *g,
            double *factr, double *pgtol, double *wa, int *iwa, char *task, int *iprint,
            char *csave, int *lsave, int *isave, double *dsave);
}

enum Activity {
    LOWER_BOUND = 0,
    UPPER_BOUND
};

class FilterAugmentedLagrangian {
    public:
        double penalty_parameter;
        double penalty_update_factor;
        double epsilon = 1e-5;
        int number_outer_iterations = 2000;
        int number_inner_iterations = 1000;
        int min_number_inner_iterations = 3;
        //double sigma = 0.01;
        bool use_eqp_step = true;
    
        FilterAugmentedLagrangian();
        
        double compute_eta(std::vector<double>& constraints);
        double compute_omega(Problem& problem, std::vector<double>& x, std::vector<double>& al_gradient);
        std::vector<double> update_multipliers(Problem& problem, std::vector<double>& c_j, std::vector<double>& constraint_multipliers);
        std::vector<double> run_restoration_phase(Problem& problem, std::vector<double>& x_jplus1, std::vector<double>& y_k, double& restoration_objective, std::vector<double>& l, std::vector<double>& u, std::vector<int>& nbd, std::vector<double>& wa, std::vector<int>& iwa);
        SubproblemSolution compute_eqp_step(Problem& problem, std::vector<double>& x_kplus1, std::vector<double>& y_kplus1, std::vector<double>& c_kplus1, std::vector<double>& l, std::vector<double>& u);
        int solve(std::string problem_name);
        void reset_bfgs();
        
        int number_evaluations_objective = 0;
        
        // initialize the AL filter
        double beta = 0.999;
        double gamma = 0.001;
        FilterConstants filter_constants = {beta, gamma};
        Filter filter;
    
    private:
        int limited_memory_size;
        /* Fortran parameters needed by lbfgsb.f */
        char task_[60];
        char csave_[60];
        int lsave_[4];
        int isave_[44];
        double dsave_[29];
        int iprint_ = -1;
        double factr_ = 1e3;
        //double pgtol_ = 1e-5;
        double pgtol_ = 1e-5;
};

FilterAugmentedLagrangian::FilterAugmentedLagrangian(): penalty_parameter(10.), penalty_update_factor(10.), filter(filter_constants), limited_memory_size(5) {
}

void FilterAugmentedLagrangian::reset_bfgs() {
    strcpy(this->task_, "START");
    std::cout << "Resetting BFGS\n";
}

// evaluate the augmented Lagrangian at x,y, where y=constraint_multipliers
double compute_al(Problem& problem, std::vector<double>& x, std::vector<double>& constraints, std::vector<double>& constraint_multipliers, double penalty_parameter) {
    // contribution of the objective
    double al = problem.objective(x);
    
    // contribution of the constraints
    for (int j = 0; j < problem.number_constraints; j++) {
        al -= constraint_multipliers[j]*constraints[j]; // f = f - lambda[j]*c[j]
        al += penalty_parameter*constraints[j]*constraints[j]/2.; // f = f + rho/2 c[j]^2
    }
    return al;
}

// evaluate the gradient of the augmented Lagrangian at x,y, where y=constraint_multipliers
std::vector<double> compute_al_gradient(Problem& problem, std::vector<double>& x, std::vector<double>& c, std::vector<double>& y, double penalty_parameter) {
    // gradient of the objective
    std::vector<double> al_gradient = problem.objective_dense_gradient(x);
    
    // Jacobian of the constraints wrt original variables
    std::vector<std::map<int,double> > constraints_sparse_jacobian = problem.constraints_sparse_jacobian(x);
    for (int j = 0; j < problem.number_constraints; j++) {
        double factor_j = y[j] - penalty_parameter*c[j];
        for (std::pair<const int, double> element: constraints_sparse_jacobian[j]) {
            int i = element.first;
            double derivative_i = element.second;
            al_gradient[i] -= factor_j*derivative_i;
        }
    }
    // gradient of the inequality constraints wrt slacks
    for (std::pair<const int, int> element: problem.inequality_constraints) {
        int j = element.first; // index of the constraint
        double derivative = y[j] - penalty_parameter*c[j];
        al_gradient.push_back(derivative);
    }
    return al_gradient;
}

// constraint violation (infeasibility)
double FilterAugmentedLagrangian::compute_eta(std::vector<double>& constraints) {
    return norm_2(constraints);
}

// residual of first-order conditions
double FilterAugmentedLagrangian::compute_omega(Problem& problem, std::vector<double>& x, std::vector<double>& al_gradient) {
    // primal variables
    double residual = 0.;
    for (int i = 0; i < problem.number_variables; i++) {
        double error = std::min(x[i] - problem.variables_bounds[i].lb, std::max(x[i] - problem.variables_bounds[i].ub, al_gradient[i]));
        residual += error*error;
    }
    // slacks
    for (std::pair<const int, int> element: problem.inequality_constraints) {
        int j = element.first; // index of the constraint
        int slack_index = problem.number_variables + element.second;
        double slack_value = x[slack_index];
        double error = std::min(slack_value - problem.constraints_bounds[j].lb, std::max(slack_value - problem.constraints_bounds[j].ub, al_gradient[slack_index]));
        residual += error*error;
    }
    return std::sqrt(residual);
}

std::vector<double> compute_constraints(Problem& problem, std::vector<double>& x) {
    std::vector<double> original_constraints = problem.evaluate_constraints(x);
    std::vector<double> constraints(problem.number_constraints);
    for (int j = 0; j < problem.number_constraints; j++) {
        try {
            // inequality constraint: subtract slack value
            int slack_index = problem.number_variables + problem.inequality_constraints.at(j);
            constraints[j] = original_constraints[j] - x[slack_index];
        }
        catch (std::out_of_range) {
            // equality constraint: subtract bound
            constraints[j] = original_constraints[j] - problem.constraints_bounds[j].lb;
        }
    }
    return constraints;
}

std::vector<double> FilterAugmentedLagrangian::update_multipliers(Problem& problem, std::vector<double>& c_j, std::vector<double>& constraint_multipliers) {
    /* compute the new multipliers by using first-order update formula: y_trial = y - rho*c */
    std::vector<double> trial_constraint_multipliers(problem.number_constraints);
    for (int j = 0; j < problem.number_constraints; j++) {
        trial_constraint_multipliers[j] = constraint_multipliers[j] - this->penalty_parameter*c_j[j];
    }
    return trial_constraint_multipliers;
}

std::vector<double> FilterAugmentedLagrangian::run_restoration_phase(Problem& problem, std::vector<double>& x_jplus1, std::vector<double>& y_k, double& restoration_objective, std::vector<double>& l, std::vector<double>& u, std::vector<int>& nbd, std::vector<double>& wa, std::vector<int>& iwa) {
    std::vector<double> x_restoration(x_jplus1);
    int n = x_restoration.size();
    std::vector<double> restoration_objective_gradient(x_restoration.size());
    
    int bfgs_iteration = 0;
    bool stop_bfgs_restoration = false;
    this->reset_bfgs();
    
    // BFGS loop (lbfgsb.f uses reverse communication to get function and gradient values)
    while (!stop_bfgs_restoration) {
        /* call L-BFGS-B */
        setulb_(&n, &this->limited_memory_size, x_restoration.data(), l.data(), u.data(), nbd.data(), &restoration_objective, restoration_objective_gradient.data(), &this->factr_, &this->pgtol_, wa.data(), iwa.data(), this->task_, &this->iprint_, this->csave_, this->lsave_, this->isave_, this->dsave_);
        
        // evaluate Augmented Lagrangian and its gradient
        if (strncmp(this->task_, "FG", 2) == 0) {
            // compute constraint violation at the point x_restoration
            std::vector<double> constraints = compute_constraints(problem, x_restoration);
            restoration_objective = 0.;
            for (int j = 0; j < problem.number_constraints; j++) {
                restoration_objective += constraints[j]*constraints[j];
            }
            // compute constraint violation gradient at the point x_restoration
            restoration_objective_gradient = std::vector<double>(problem.number_variables);
            for (int j = 0; j < problem.number_constraints; j++) {
                double factor_j = 2.*constraints[j];
                std::vector<double> constraint_j_gradient = problem.constraint_dense_gradient(j, x_restoration);
                for (int i = 0; i < problem.number_variables; i++) {
                    restoration_objective_gradient[i] += factor_j*constraint_j_gradient[i];
                }
            }
            // gradient of the inequality constraints wrt slacks
            for (std::pair<const int, int> element: problem.inequality_constraints) {
                int j = element.first; // index of the constraint
                restoration_objective_gradient.push_back(-2.*constraints[j]);
            }
        }
        // new x value
        else if (strncmp(this->task_, "NEW_X", 5) == 0) {
            std::cout << "Restoration iteration " << bfgs_iteration << "\n";
            std::cout << "x = "; print_vector(std::cout, x_restoration);
            std::cout << "f = " << restoration_objective << "\n";
            std::cout << "g = "; print_vector(std::cout, restoration_objective_gradient);
            bfgs_iteration++;
            
            std::vector<double> c = compute_constraints(problem, x_restoration);
            std::vector<double> al_gradient = compute_al_gradient(problem, x_restoration, c, y_k, this->penalty_parameter);
            double omega = this->compute_omega(problem, x_restoration, al_gradient);
            double eta = this->compute_eta(c);
            
            if (filter.accept(eta, omega)) {
                std::cout << "Filter-acceptable\n";
                stop_bfgs_restoration = true;
            }
        }
        // termination criteria
        else if (strncmp(this->task_, "CONV", 4) == 0) {
            std::cout << "BFGS converged with stationary point\n";
            stop_bfgs_restoration = true;
        }
        else if (this->number_inner_iterations < bfgs_iteration) {
            throw std::logic_error("BFGS reached max iter in restoration");
        }
        else if (strncmp(this->task_, "ABNO", 4) == 0 || strncmp(this->task_, "ERROR", 5) == 0) {
            std::cout << "Task: " << this->task_ << "\n";
            throw std::logic_error("Error in BFGS");
        } 
    }
    this->reset_bfgs();
    std::cout << "BFGS restoration phase exited with solution x*: "; print_vector(std::cout, x_restoration);
    return x_restoration;
}

// constraints contains the evaluation of the constraints of the slacked problem
SubproblemSolution FilterAugmentedLagrangian::compute_eqp_step(Problem& problem, std::vector<double>& x_kplus1, std::vector<double>& y_kplus1, std::vector<double>& c_kplus1, std::vector<double>& l, std::vector<double>& u) {
    std::cout << "************** SOC **************\n";

    // compute the Hessian of the Lagrangian
    std::vector<double> y_al(problem.number_constraints);
    for (int j = 0; j < problem.number_constraints; j++) {
        y_al[j] = y_kplus1[j] - this->penalty_parameter*c_kplus1[j];
    }
    CSCMatrix hessian = problem.lagrangian_hessian(x_kplus1, problem.objective_sign, y_al);
    
    /* constraint Jacobian */
    std::vector<std::map<int, double> > constraints_jacobian = problem.constraints_sparse_jacobian(x_kplus1);
    // add the contribution of the slacks
    for (std::pair<const int, int> constraint: problem.inequality_constraints) {
        int j = constraint.first;
        int slack_index = problem.number_variables + constraint.second;
        constraints_jacobian[j][slack_index] = -1.;
    }

    // Hessian of the AL
    ArgonotMatrix al_hessian = hessian.to_ArgonotMatrix(x_kplus1.size());
    
    /* perform matrix addition: Hessian(x, y - rho*c) + rho \nabla c \nabla c^T */
    for (int j = 0; j < problem.number_constraints; j++) {
        // add to KKT_matrix all the outer products
        for (std::pair<const int, double> row_element: constraints_jacobian[j]) {
            int row_index = row_element.first;
            double row_derivative = row_element.second;
            for (std::pair<const int, double> column_element: constraints_jacobian[j]) {
                int column_index = column_element.first;
                double column_derivative = column_element.second;
                // upper triangular matrix
                if (column_index >= row_index) {
                    // add product of components
                    al_hessian.add_term(this->penalty_parameter*row_derivative*column_derivative, row_index, column_index);
                }
            }
        }
    }
    // convert to CSC matrix
    CSCMatrix csc_hessian = al_hessian.to_CSC();
    
    // fix bounds of variables in active set
    std::vector<Range> variables_bounds(x_kplus1.size());
    for (unsigned int i = 0; i < x_kplus1.size(); i++) {
        if (x_kplus1[i] == l[i]) {
            std::cout << "Active x" << i << "\n";
            variables_bounds[i] = {l[i] - x_kplus1[i], l[i] - x_kplus1[i]};
        }
        else if (x_kplus1[i] == u[i]) {
            std::cout << "Active x" << i << "\n";
            variables_bounds[i] = {u[i] - x_kplus1[i], u[i] - x_kplus1[i]};
        }
        else {
            variables_bounds[i] = {l[i] - x_kplus1[i], u[i] - x_kplus1[i]};
        }
    }
    //print_vector(std::cout, );

    // objective gradient
    std::vector<double> al_gradient = compute_al_gradient(problem, x_kplus1, c_kplus1, y_al, this->penalty_parameter);
    std::map<int, double> linear_objective;
    for (unsigned int i = 0; i < al_gradient.size(); i++) {
        linear_objective[i] = al_gradient[i];
    }

    // set up bounds of linearized constraints
    std::vector<Range> constraints_bounds(problem.number_constraints);
    for (int j = 0; j < problem.number_constraints; j++) {
        constraints_bounds[j] = {-c_kplus1[j], -c_kplus1[j]};
    }
    
    // solve the QP
    BQPDSolver solver(csc_hessian.column_start, csc_hessian.row_number);
    solver.allocate(x_kplus1.size(), problem.number_constraints);
    std::vector<double> d0(x_kplus1.size());
    SubproblemSolution eqp_step = solver.solve_QP(variables_bounds, constraints_bounds, linear_objective, constraints_jacobian, csc_hessian, d0);
    
    // compute multiplier displacement
    for (int j = 0; j < problem.number_constraints; j++) {
        eqp_step.multipliers.constraints[j] -= y_kplus1[j];
    }
    
    return eqp_step;
}

int FilterAugmentedLagrangian::solve(std::string problem_name) {
    // create the problem
    AMPLModel problem = AMPLModel(problem_name);

    // initial primal and dual points
    std::vector<double> x_k = problem.primal_initial_solution();
    std::vector<double> y_k = problem.dual_initial_solution();
    
    // add slacks as primal variables
    for (std::pair<const int, int> element: problem.inequality_constraints) {
        int j = element.first; // index of the constraint
        x_k.push_back(problem.evaluate_constraint(j, x_k));
    }
    
    // compute bounds and variable status
    int n = x_k.size();
    std::vector<double> l(n);
    std::vector<double> u(n);
    std::vector<ConstraintType> variable_status(n);
    for (int i = 0; i < problem.number_variables; i++) {
        l[i] = problem.variables_bounds[i].lb;
        u[i] = problem.variables_bounds[i].ub;
        variable_status[i] = problem.variable_status[i];
    }
    for (std::pair<const int, int> element: problem.inequality_constraints) {
        int j = element.first;
        int slack_index = problem.number_variables + element.second;
        l[slack_index] = problem.constraints_bounds[j].lb;
        u[slack_index] = problem.constraints_bounds[j].ub;
        variable_status[slack_index] = problem.constraint_status[j];
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
    
    // evaluate the initial point
    std::vector<double> c_k = compute_constraints(problem, x_k);
    double al_k = compute_al(problem, x_k, c_k, y_k, this->penalty_parameter);
    std::vector<double> al_gradient_k = compute_al_gradient(problem, x_k, c_k, y_k, this->penalty_parameter);
    double omega_k = this->compute_omega(problem, x_k, al_gradient_k);
    double eta_k = this->compute_eta(c_k);
    this->number_evaluations_objective++;
    
    std::cout << "Initial x: "; print_vector(std::cout, x_k);
    std::cout << "Initial constraints: "; print_vector(std::cout, c_k);
    std::cout << "Initial measures: " << eta_k << " " << omega_k << "\n";
    
    // initialize the filter with initial point's entries
    if (0. < eta_k) {
        filter.add(eta_k, omega_k);
        filter.upper_bound = std::max(100., 1.25*eta_k);
    }
    
    //filter.upper_bound = delta*std::max(1., 1.25*eta_k);
    
    /* memory allocation for L-BFGS-B (limited memory & factors allocation) */
    std::vector<double> wa(this->limited_memory_size*(2*n + 11*this->limited_memory_size + 8) + 5*n);
    std::vector<int> iwa(3*n);

    // optimization loop
    bool convergence = false;
    int iterations = 0;
    std::string termination_message;
    this->reset_bfgs();
    
    while (!convergence && iterations < number_outer_iterations) {
        bool filter_acceptable = false;
        std::vector<double> x_j(x_k);
        std::vector<double> y_j(y_k);
        std::vector<double> c_j(c_k);
        double al_j = al_k;
        double omega_j = omega_k;
        double eta_j = eta_k;
        
        std::vector<double> x_jplus1(x_j);
        std::vector<double> y_jplus1(y_j);
        std::vector<double> c_jplus1(c_j);
        double al_jplus1 = al_j;
        std::vector<double> al_gradient_jplus1(al_gradient_k);
        double omega_jplus1 = omega_j;
        double eta_jplus1 = eta_j;
        
        while (!filter_acceptable && !convergence) {
            if (strncmp(this->task_, "START", 5) == 0) {
                std::cout << "\n## Starting BFGS from scratch from ";
            }
            else {
                std::cout << "\n## Warm-starting BFGS from ";
            }
            print_vector(std::cout, x_j);
            std::cout << "Filter upper bound is " << filter.upper_bound << "\n";
            std::cout << "Penalty parameter is " << this->penalty_parameter << "\n";
            std::cout << "Fixed multipliers: "; print_vector(std::cout, y_k);
            
            /************************ start BFGS**********************/
            // approximately minimize Augmented Lagrangian subproblem
            bool sufficient_progress = false;
            int bfgs_iteration = 0;
            
            // BFGS loop (reverse communication) to find a point that guarantees sufficient decrease
            while (!sufficient_progress) {
                /* call L-BFGS-B */
                setulb_(&n, &this->limited_memory_size, x_jplus1.data(), l.data(), u.data(), nbd.data(), &al_jplus1, al_gradient_jplus1.data(), &this->factr_, &this->pgtol_, wa.data(), iwa.data(), this->task_, &this->iprint_, this->csave_, this->lsave_, this->isave_, this->dsave_);
                
                // convergence to a stationary point
                if (strncmp(this->task_, "CONV", 4) == 0) {
                    std::cout << "BFGS converged with stationary point\n";
                    c_jplus1 = compute_constraints(problem, x_jplus1);
                    al_jplus1 = compute_al(problem, x_jplus1, c_jplus1, y_k, this->penalty_parameter);
                    al_gradient_jplus1 = compute_al_gradient(problem, x_jplus1, c_jplus1, y_k, this->penalty_parameter);
                    omega_jplus1 = this->compute_omega(problem, x_jplus1, al_gradient_jplus1);
                    eta_jplus1 = this->compute_eta(c_jplus1);
                    y_jplus1 = this->update_multipliers(problem, c_jplus1, y_k);
                    // objective was already computed, do not count evaluation
                    
                    sufficient_progress = true;
                }
                // possible error cases
                else if (strncmp(this->task_, "ERROR", 5) == 0) {
                    std::cout << this->task_ << "\n";
                    throw std::logic_error("Error in BFGS");
                }
                else if (this->number_inner_iterations <= bfgs_iteration) {
                    throw std::logic_error("BFGS reached max iter");
                }
                else if (strncmp(this->task_, "ABNO", 4) == 0) {
                    std::cout << this->task_ << "\n";
                    throw std::logic_error("BFGS: abnormality in line-search");
                }
                else if (strncmp(this->task_, "FG", 2) == 0) {
                    // compute augmented Lagrangian and its gradient at the point (x_trial, y)
                    c_jplus1 = compute_constraints(problem, x_jplus1);
                    al_jplus1 = compute_al(problem, x_jplus1, c_jplus1, y_k, this->penalty_parameter);
                    al_gradient_jplus1 = compute_al_gradient(problem, x_jplus1, c_jplus1, y_k, this->penalty_parameter);
                    omega_jplus1 = this->compute_omega(problem, x_jplus1, al_gradient_jplus1);
                    eta_jplus1 = this->compute_eta(c_jplus1);
                    y_jplus1 = this->update_multipliers(problem, c_jplus1, y_k);
                    this->number_evaluations_objective++;
                }
                else if (strncmp(this->task_, "NEW_X", 5) == 0) {
                    std::cout << "Iteration " << bfgs_iteration << "\n";
                    
                    c_jplus1 = compute_constraints(problem, x_jplus1);
                    omega_jplus1 = this->compute_omega(problem, x_jplus1, al_gradient_jplus1);
                    eta_jplus1 = this->compute_eta(c_jplus1);
                    y_jplus1 = this->update_multipliers(problem, c_jplus1, y_k);
                    // objective was already computed, do not count evaluation
                    
                    std::cout << "  x = "; print_vector(std::cout, x_jplus1);
                    std::cout << "  c(x) = "; print_vector(std::cout, c_jplus1);
                    std::cout << "  y = "; print_vector(std::cout, y_jplus1);
                    std::cout << "  f = " << al_jplus1 << "\n";
                    std::cout << "  g = "; print_vector(std::cout, al_gradient_jplus1);
                    std::cout << "  η = " << eta_jplus1 << ", ω = " << omega_jplus1 << "\n";
                        
                    bfgs_iteration++;
                    // perform at least min_number_inner_iterations iterations before testing convergence
                    if (min_number_inner_iterations <= bfgs_iteration) {
                        if (filter.accept(eta_jplus1, omega_jplus1)) {
                            std::cout << "Filter-acceptable\n";
                            sufficient_progress = true;
                        }
                    }
                }
            }
            /************************ end BFGS**********************/
            // we have here sufficient decrease
            
            std::cout << "End of BFGS loop\n";
            std::cout << "Measures: (η = " << eta_jplus1 << ", ω = " << omega_jplus1 << ")\n";
            std::cout << "There are " << filter.entries.size() << " entries in the filter\n";
            
            bool switch_to_restoration = false;
            // Eq 3.15
            if (filter.entries.size() > 0 && eta_jplus1 >= beta*std::max(filter.omega_min()/gamma, filter.eta_min())) {
                switch_to_restoration = true;
            }
            // Eq 3.16
            else if (filter.entries.size() > 0 && omega_jplus1 <= epsilon && eta_jplus1 >= beta*filter.eta_min()) { 
                switch_to_restoration = true;
            }
            else if (strncmp(this->task_, "CONV", 4) == 0 && !filter.accept(eta_j, omega_j)) {
                // convergence to a stationary point, but no progress made towards optimality
                switch_to_restoration = true;
            }
            
            // restoration phase
            if (switch_to_restoration) { 
                std::cout << "  *** Switching to restoration phase to find filter-acceptable point\n";

                // parameter update
                std::cout << "Increasing parameter\n";
                this->penalty_parameter *= this->penalty_update_factor;
                
                // run restoration phase
                double restoration_objective;
                x_jplus1 = run_restoration_phase(problem, x_jplus1, y_k, restoration_objective, l, u, nbd, wa, iwa);
                
                // if convergence towards a stationary point
                if (restoration_objective < epsilon) {
                    c_jplus1 = compute_constraints(problem, x_jplus1);
                    al_jplus1 = compute_al(problem, x_jplus1, c_jplus1, y_k, this->penalty_parameter);
                    al_gradient_jplus1 = compute_al_gradient(problem, x_jplus1, c_jplus1, y_k, this->penalty_parameter);
                    omega_jplus1 = this->compute_omega(problem, x_jplus1, al_gradient_jplus1);
                    eta_jplus1 = this->compute_eta(c_jplus1);
                    y_jplus1 = this->update_multipliers(problem, c_jplus1, y_k);
                    this->number_evaluations_objective++;

                    std::cout << "Restoration c(x*) = "; print_vector(std::cout, c_jplus1);
                    std::cout << "Restoration measures: (η = " << eta_jplus1 << ", ω = " << omega_jplus1 << ")\n";
                }
                else {
                    termination_message = "Restoration phase ended with non-feasible FJ point";
                    convergence = true;
                }
            }
            else {
                // carry on with the optimization
            }
            
            // update j <- j+1
            x_j = x_jplus1;
            y_j = y_jplus1;
            c_j = c_jplus1;
            al_j = al_jplus1;
            eta_j = eta_jplus1;
            omega_j = omega_jplus1;
            
            // test filter acceptance
            filter_acceptable = filter.accept(eta_j, omega_j);
            if (filter_acceptable) {
                std::cout << "The point is filter-acceptable\n";
            }
            else {
                std::cout << "The point is NOT filter-acceptable\n";
            }
        }
        std::cout << "Exiting inner loop\n";
        // here, we have a filter-acceptable point
        std::vector<double> x_kplus1(x_j);
        std::vector<double> y_kplus1(y_j);
        std::vector<double> c_kplus1(c_j);
        double al_kplus1(al_j);
        double eta_kplus1(eta_j);
        double omega_kplus1(omega_j);
        iterations++;
        // since the multipliers change, we reset BFGS
        this->reset_bfgs();
        
        if (eta_j <= epsilon && omega_j <= epsilon) {
            convergence = true;
        }
        
        /* figure out active set */
        bool inactive_variables = false;
        for (unsigned int i = 0; i < x_kplus1.size(); i++) {
            if (l[i] < x_kplus1[i] || x_kplus1[i] < u[i]) {
                inactive_variables = true;
                break;
            }
        }
        if (this->use_eqp_step && !convergence && inactive_variables) {
            // optional: perform second-order correction + LS + recompute eta and omega
            SubproblemSolution eqp_step = compute_eqp_step(problem, x_kplus1, y_kplus1, c_kplus1, l, u);
            std::cout << "EQP step:\n";
            std::cout << "delta x = "; print_vector(std::cout, eqp_step.x);
            std::cout << "delta constraint_multipliers = "; print_vector(std::cout, eqp_step.multipliers.constraints);

            // run line search along EQP step
            double step_length = 1.;
            double min_step_length = 0.01;
            bool is_accepted = false;
            bool max_iterations_reached = false;
            while(!is_accepted && !max_iterations_reached) {
                try {
                    // take a step of length step_length along EQP step
                    std::vector<double> x_j_soc = add_vectors(x_kplus1, eqp_step.x, step_length);
                    std::vector<double> y_j_soc = add_vectors(y_kplus1, eqp_step.multipliers.constraints, step_length);

                    // compute measures
                    std::vector<double> c_j_soc = compute_constraints(problem, x_j_soc);
                    double al_j_soc = compute_al(problem, x_j_soc, c_j_soc, y_j_soc, this->penalty_parameter);
                    std::vector<double> al_gradient_j_soc = compute_al_gradient(problem, x_j_soc, c_j_soc, y_k, this->penalty_parameter);
                    double eta_j_soc = compute_eta(c_j_soc);
                    double omega_j_soc = compute_omega(problem, x_j_soc, al_gradient_j_soc);
                    this->number_evaluations_objective++;
                    
                    // accept the point if filter-acceptable
                    //  && filter.improves_current_iterate(eta_j, omega_j, eta_j_soc, omega_j_soc)
                    if (filter.accept(eta_j_soc, omega_j_soc)) {
                        is_accepted = true;
                        std::cout << "Accepting point with step length " << step_length << "\n";
                        std::cout << "Current measures: (η = " << eta_j << ", ω = " << omega_j << ")\n";
                        std::cout << "SOC measures: (η = " << eta_j_soc << ", ω = " << omega_j_soc << ")\n";
                        // update the primal-dual point with the SOC solution
                        // update k <- k+1
                        x_kplus1 = x_j_soc;
                        y_kplus1 = y_j_soc;
                        c_kplus1 = c_j_soc;
                        al_kplus1 = al_j_soc;
                        eta_kplus1 = eta_j_soc;
                        omega_kplus1 = omega_j_soc;
                    }
                    else if (step_length < min_step_length) {
                        std::cout << "Rejecting point in line-search\n";
                        max_iterations_reached = true;
                        // keep the filter-acceptable (x_kplus1, y_kplus1)
                    }
                }
                catch (const std::invalid_argument& e) {
                    std::cout << "IEEE error while evaluating. Backtracking\n";
                    is_accepted = false;
                }
                
                if (!is_accepted && !max_iterations_reached) {
                    step_length /= 2.;
                }
            }
            // end of SOC line search
        }
        else {
            x_kplus1 = x_j;
            y_kplus1 = y_j;
            c_kplus1 = c_j;
            al_kplus1 = al_j;
            eta_kplus1 = eta_j;
            omega_kplus1 = omega_j;
        }
        
        if (0. < eta_kplus1) {
            std::cout << "Adding entries to the filter\n";
            filter.add(eta_kplus1, omega_kplus1);
        }
        x_k = x_kplus1;
        y_k = y_kplus1;
        c_k = c_kplus1;
        al_k = al_kplus1;
        eta_k = eta_kplus1;
        omega_k = omega_kplus1;
        
        if (eta_k <= epsilon && omega_k <= epsilon) {
            convergence = true;
        }
    }
    if (convergence) {
        termination_message = "Subsolver converged within given tolerance";
    }
    else if (iterations >= number_outer_iterations) {
        termination_message = "number_outer_iterations was reached";
    }

    std::cout << "\nOptimization summary:\n";
    std::cout << termination_message << "\n";
    std::cout << "x = "; print_vector(std::cout, x_k);
    std::cout << "constraint multipliers = "; print_vector(std::cout, y_k);
    std::cout << "Objective evaluations: " << this->number_evaluations_objective << "\n";
    std::cout << "Constraint evaluations: " << problem.number_eval_constraints << "\n";
    std::cout << "Outer iterations: " << iterations << "\n";
    
    return 0;
}

int main(int argc, char* argv[]) {
    std::string problem_name = std::string(argv[argc - 1]);
    FilterAugmentedLagrangian filter_al;
    return filter_al.solve(problem_name);
}
