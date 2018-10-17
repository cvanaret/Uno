#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <string.h>
#include "LBFGSB.hpp"
#include "Utils.hpp"

//extern "C" {
//    #include "lbfgsb.h"
//}

extern "C" {

    void setulb_(int *n, int *m, double *x, double *l, double *u, int *nbd, double *f, double *g,
            double *factr, double *pgtol, double *wa, int *iwa, char *task, int *iprint,
            char *csave, int *lsave, int *isave, double *dsave);//, long int, long int);
}

LBFGSB::LBFGSB(int memory_size): rho(200.), memory_size(memory_size) {
}

void LBFGSB::initialize(std::map<int,int> slacked_constraints) {
    this->slacked_constraints_ = slacked_constraints;
    return;
}

LocalSolution LBFGSB::solve(Problem& problem, Iterate& current_iterate) {
    std::vector<double> x(current_iterate.x);
    /* determine the bounds of the variables + slacks */
    int n = problem.number_variables + this->slacked_constraints_.size();
    std::vector<int> nbd(n);
    std::vector<double> l(n);
    std::vector<double> u(n);
    for (int i = 0; i < problem.number_variables; i++) {
        l[i] = problem.variable_lb[i];
        u[i] = problem.variable_ub[i];
        if (problem.variable_status[i] == UNBOUNDED) {
            nbd[i] = 0;
        }
        else if (problem.variable_status[i] == BOUNDED_LOWER) {
            nbd[i] = 1;
        }
        else if (problem.variable_status[i] == BOUNDED_UPPER) {
            nbd[i] = 3;
        }
        else {
            nbd[i] = 2;
        }
    }
    for (std::pair<const int, int> element: this->slacked_constraints_) {
        int j = element.first;
        int current_slack = element.second;
        l[problem.number_variables + current_slack] = problem.constraint_lb[j];
        u[problem.number_variables + current_slack] = problem.constraint_ub[j];
        
        if (problem.constraint_status[j] == UNBOUNDED) {
            nbd[problem.number_variables + current_slack] = 0;
        }
        else if (problem.constraint_status[j] == BOUNDED_LOWER) {
            nbd[problem.number_variables + current_slack] = 1;
        }
        else if (problem.constraint_status[j] == BOUNDED_UPPER) {
            nbd[problem.number_variables + current_slack] = 3;
        }
        else {
            nbd[problem.number_variables + current_slack] = 2;
        }
    }
    
    /* memory allocation */
    std::vector<double> wa(this->memory_size*(2*n + 11*this->memory_size + 8) + 5*n);
    std::vector<int> iwa(3*n);
    
    // optimization loop
    double f;
    std::vector<double> g(n);
    
    strcpy(this->task_, "START");
    while (strncmp(this->task_, "FG", 2) == 0 || strncmp(this->task_, "NEW_X", 5) == 0 || strncmp(this->task_, "START", 5) == 0) {
        /* call L-BFGS-B */
        setulb_(&n, &this->memory_size, x.data(), l.data(), u.data(), nbd.data(), &f, g.data(), &this->factr_, &this->pgtol_, wa.data(), iwa.data(),
                this->task_, &this->iprint_, this->csave_, this->lsave_, this->isave_, this->dsave_);//, (long int) 60, (long int) 60);
        
        std::cout << "Current task: " << this->task_ << "\n";
        
        // evaluate Augmented Lagrangian and its gradient
        if (strncmp(this->task_, "FG", 2) == 0) {
            std::cout << "x: "; print_vector(std::cout, x);
            std::vector<double> constraints = problem.evaluate_constraints(x);
            f = this->compute_augmented_lagrangian_(problem, x, constraints, current_iterate.constraint_multipliers);
            g = this->compute_augmented_lagrangian_gradient_(problem, x, constraints, current_iterate.constraint_multipliers);
            std::cout << "f is " << f << "\n";
            std::cout << "g is "; print_vector(std::cout, g);
        }
    }

    //! ... print output of this example
    std::cout << "Final L-BFGS-B Solution\n";
    std::cout << "lower bound   x-value      upper bound  gradient\n";
    double reduced_gradient = 0.;
    for (int i = 0; i < n; i++) {
        std::cout << x[i] << " in [" << l[i] << ", " << u[i] << "]\tderivative: " << g[i] << "\n";
        reduced_gradient += std::abs(std::min(x[i] - l[i], u[i] - x[i]) * g[i]);
    }; // end for
    std::cout << "Reduced Gradient Norm = " << reduced_gradient << "\n";
    
    /* compute the new multipliers */
    std::vector<double> constraints = problem.evaluate_constraints(x);
    std::vector<double> constraint_multipliers(current_iterate.constraint_multipliers);
    for (int j = 0; j < problem.number_constraints; j++) {
        double constraint_value;
        try {
            // inequality constraint
            int current_slack = this->slacked_constraints_[j];
            constraint_value = constraints[j] - x[problem.number_variables + current_slack];
        }
        catch (std::out_of_range) {
            // equality constraint
            constraint_value = constraints[j] - problem.constraint_lb[j];
        }
        constraint_multipliers[j] -= this->rho*constraint_value;
    }
    LocalSolution solution(x, g, constraint_multipliers);
    solution.status = OPTIMAL;

    return solution;
}

double LBFGSB::compute_augmented_lagrangian_(Problem& problem, std::vector<double>& x, std::vector<double>& constraints, std::vector<double>& constraint_multipliers) {
    // contribution of the objective
    double f = problem.objective(x);
    // contribution of the constraints
    for (int j = 0; j < problem.number_constraints; j++) {
        double constraint_value;
        try {
            // inequality constraint
            int current_slack = this->slacked_constraints_[j];
            constraint_value = constraints[j] - x[problem.number_variables + current_slack];
        }
        catch (std::out_of_range) {
            // equality constraint
            constraint_value = constraints[j] - problem.constraint_lb[j];
        }
        f -= constraint_multipliers[j]*constraint_value;
        f += this->rho/2.*constraint_value*constraint_value;
    }
    return f;
}

std::vector<double> LBFGSB::compute_augmented_lagrangian_gradient_(Problem& problem, std::vector<double>& x, std::vector<double>& constraints, std::vector<double>& constraint_multipliers) {
    // gradient of the objective
    std::vector<double> augmented_lagrangian_gradient = problem.objective_dense_gradient(x);
    // gradient of the constraints wrt the variablescurrent_iterate
    for (int j = 0; j < problem.number_constraints; j++) {
        double constraint_value;
        try {
            // inequality constraint
            int current_slack = this->slacked_constraints_[j];
            constraint_value = constraints[j] - x[problem.number_variables + current_slack];
        }
        catch (std::out_of_range) {
            // equality constraint
            constraint_value = constraints[j] - problem.constraint_lb[j];
        }
        double factor = constraint_multipliers[j] - this->rho*constraint_value;
        // update the gradient
        std::vector<double> constraint_gradient = problem.constraint_dense_gradient(j, x);
        for (int i = 0; i < problem.number_variables; i++) {
            augmented_lagrangian_gradient[i] -= factor*constraint_gradient[i];
        }
    }
    // gradient of the constraints wrt the slacks
    for (std::pair<const int, int> element: slacked_constraints_) {
        int j = element.first;
        int current_slack = element.second;
        double constraint_value = constraints[j] - x[problem.number_variables + current_slack];
        double derivative = constraint_multipliers[j] - this->rho*constraint_value;
        augmented_lagrangian_gradient.push_back(derivative);
    }
    return augmented_lagrangian_gradient;
}