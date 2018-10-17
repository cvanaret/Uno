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

LBFGSB::LBFGSB(): rho(10.) {
}

void LBFGSB::initialize(std::map<int,int> slacked_constraints) {
    this->slacked_constraints = slacked_constraints;
    return;
}

LocalSolution LBFGSB::solve(Problem& problem, Iterate& current_iterate) {
    std::vector<double> x(current_iterate.x);
    // determine the bounds of the variables + slacks
    std::cout << "PROBLEM HAS " << problem.number_variables << " vars\n";
    int n = problem.number_variables + this->slacked_constraints.size();
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
    for (std::pair<const int, int> element: this->slacked_constraints) {
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
    
    int m = 5, iprint = 1;
    double factr = 1e5, pgtol = 1e-5;
    char task[60];
    char csave[60];
    int lsave[4];
    int isave[44];
    double dsave[29];

    std::vector<double> wa(2*m*n + 11*m*m + 5*n + 8*m);
    std::vector<int> iwa(3 * n);
    
    // optimization loop
    double f;
    std::vector<double> g(n);
    strcpy(task, "START");
    while (strncmp(task, "FG", 2) == 0 || strncmp(task, "NEW_X", 5) == 0 || strncmp(task, "START", 5) == 0) {
        // call L-BFGS-B
        setulb_(&n, &m, x.data(), l.data(), u.data(), nbd.data(), &f, g.data(), &factr, &pgtol, wa.data(), iwa.data(),
                task, &iprint, csave, lsave, isave, dsave);//, (long int) 60, (long int) 60);
        
        // evaluate Augmented Lagrangian and its gradient
        if (strncmp(task, "FG", 2) == 0) {
            std::cout << "x: "; print_vector(std::cout, x);
            f = this->compute_augmented_lagrangian(problem, current_iterate, x);
            g = this->compute_augmented_lagrangian_gradient(problem, current_iterate, x);
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

    std::vector<double> constraint_multipliers;
    LocalSolution solution(x, g, constraint_multipliers);
    solution.status = OPTIMAL;

    return solution;
}

double LBFGSB::compute_augmented_lagrangian(Problem& problem, Iterate& current_iterate, std::vector<double>& x) {
    // contribution of the objective
    double f = problem.objective(x);
    // contribution of the constraints
    for (int j = 0; j < problem.number_constraints; j++) {
        double constraint_value;
        try {
            // inequality constraint
            int current_slack = this->slacked_constraints[j];
            constraint_value = problem.evaluate_constraint(j, x) - current_iterate.x[problem.number_variables + current_slack];
        }
        catch (std::out_of_range) {
            // equality constraint
            constraint_value = problem.evaluate_constraint(j, x) - problem.constraint_lb[j];
        }
        f -= current_iterate.constraint_multipliers[j]*constraint_value;
        f += this->rho/2.*constraint_value*constraint_value;
    }
    return f;
}

std::vector<double> LBFGSB::compute_augmented_lagrangian_gradient(Problem& problem, Iterate& current_iterate, std::vector<double>& x) {
    std::vector<double> constraints = problem.evaluate_constraints(x);
    
    // gradient of the objective
    std::vector<double> augmented_lagrangian_gradient = problem.objective_dense_gradient(x);
    // gradient of the constraints wrt the variables
    for (int j = 0; j < problem.number_constraints; j++) {
        double constraint_value;
        try {
            // inequality constraint
            int current_slack = this->slacked_constraints[j];
            constraint_value = constraints[j] - current_iterate.x[problem.number_variables + current_slack];
        }
        catch (std::out_of_range) {
            // equality constraint
            constraint_value = constraints[j] - problem.constraint_lb[j];
        }
        double factor = this->rho*constraint_value - current_iterate.constraint_multipliers[j];
        // update the gradient
        std::vector<double> constraint_gradient = problem.constraint_dense_gradient(j, x);
        for (int i = 0; i < problem.number_variables; i++) {
            augmented_lagrangian_gradient[i] += factor*constraint_gradient[i];
        }
    }
    // gradient of the constraints wrt the slacks
    for (std::pair<const int, int> element: slacked_constraints) {
        int j = element.first;
        int current_slack = element.second;
        double constraint_value = constraints[j] - current_iterate.x[problem.number_variables + current_slack];
        double derivative = current_iterate.constraint_multipliers[j] - this->rho*constraint_value;
        augmented_lagrangian_gradient.push_back(derivative);
    }
    return augmented_lagrangian_gradient;
}