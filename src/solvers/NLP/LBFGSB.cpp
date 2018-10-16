#include <iostream>
#include <vector>
#include <cmath>
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

LocalSolution LBFGSB::solve(Problem& problem, Iterate& current_point) {
    int n = problem.number_variables;
    int m = 5, iprint = 1;
    double factr = 1e5, pgtol = 1e-5;
    char task[60];
    char csave[60];
    int lsave[4];
    int isave[44];
    double dsave[29];

    // wa is a double precision working array of length (2*m + 5)*n + 12*m^2 + 12*m
    // iwa is an integer working array of length 3nmax.
    //std::vector<int> iwa(2*m*n + 11*m*m + 5*n + 8*m);
    std::vector<double> wa(2*m*n + 11*m*m + 5*n + 8*m);
    std::cout << "wa allocated with size " << (2*m*n + 11*m*m + 5*n + 8*m) << "\n";
    std::vector<int> iwa(3 * n);
    std::cout << "iwa allocated with size " << (3 * n) << "\n";
    
    
    std::vector<int> nbd(n);
    std::vector<double> l(problem.variable_lb);
    std::vector<double> u(problem.variable_ub);
    std::vector<double> x(current_point.x);

    for (int i = 0; i < n; i++) {
        // set nbd(i) = {0|1|2|3} = {unbdd|LbdOnly|BothBnds|UbdOnly}
        // EQUAL_BOUNDS, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES, UNBOUNDED
        if (problem.variable_status[i] == UNBOUNDED) {
            nbd[i] = 0;
        }
        else if (problem.variable_status[i] == BOUNDED_LOWER) {
            nbd[i] = 1;
        }
        else if (problem.variable_status[i] == BOUNDED_BOTH_SIDES || problem.variable_status[i] == EQUAL_BOUNDS) {
            nbd[i] = 2;
        }
        else {
            nbd[i] = 3;
        }
    }

    //! ... we start the iteration by initializing task.
    strcpy(task, "START");

    // penalty parameter
    double rho = 10.;
    
    double f;
    // gradient
    std::vector<double> g(n);
    // optimization loop
    while (strncmp(task, "FG", 2) == 0 || strncmp(task, "NEW_X", 5) == 0 || strncmp(task, "START", 5) == 0) {
        // call the L-BFGS-B code.
        setulb_(&n, &m, x.data(), l.data(), u.data(), nbd.data(), &f, g.data(),
                &factr, &pgtol, wa.data(), iwa.data(), task, &iprint, csave, lsave,
                isave, dsave);//, (long int) 60, (long int) 60);
        std::cout << "x: "; print_vector(std::cout, x);
        
        // evaluate objective and gradient
        if (strncmp(task, "FG", 2) == 0) {
            // Augmented Lagrangian
            f = problem.objective(x);
            for (int j = 0; j < problem.number_constraints; j++) {
                f -= current_point.constraint_multipliers[j]*();
            }
            
            
            g = problem.objective_dense_gradient(x);
        }
    }

    //! ... print output of this example
    std::cout << "Final L-BFGS-B Solution\n";
    std::cout << "lower bound   x-value      upper bound  gradient\n";
    double reduced_gradient = 0.;
    for (int i = 0; i < n; i++) {
        std::cout << x[i] << " in [" << l[i] << ", " << u[i] << "]\tderivative:" << g[i] << "\n";
        reduced_gradient += std::abs(std::min(x[i] - l[i], u[i] - x[i]) * g[i]);
    }; // end for
    std::cout << "Reduced Gradient Norm = " << reduced_gradient << "\n";

    std::vector<double> constraint_multipliers;
    LocalSolution solution(x, g, constraint_multipliers);

    return solution;
}