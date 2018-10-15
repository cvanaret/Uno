#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include "LBFGSB.hpp"
#include <string.h>

//extern "C" {
//    #include "lbfgsb.h"
//}

extern "C" {
    
    void setulb_(int *n, int *m, double *x, double *l, double *u, int *nbd, double *f, double *g,
        double *factr, double *pgtol, double *wa, int *iwa, char *task, int *iprint,
        char *csave, int *lsave, int *isave, double *dsave, long int, long int);
}

LocalSolution LBFGSB::solve(Problem& problem, Iterate& current_point) {
    //! ... set up storage needed by L-BFGS-B solver
    int n = 10, m = 5, iprint = 1;
    double factr = 1e5, pgtol = 1e-5;
    char task[60];
    char csave[60];
    int lsave[4];
    int isave[44];
    double dsave[29];

    // wa is a double precision working array of length (2mmax + 5)nmax + 12mmax^2 + 12mmax.
    //iwa is an integer working array of length 3nmax.
    std::vector<int> iwa((2*m + 5)*n + 12*m*m + 12*m);
    std::vector<double> wa(3*n);

    //! ... copy/define lower/upper bounds & set nbd(i) = {0|1|2|3} = {unbdd|LbdOnly|BothBnds|UbdOnly}
    std::vector<int> nbd(n);
    std::vector<double> l(n);
    std::vector<double> u(n);

    for (int i = 0; i < n; i++) {
        nbd[i] = 2;
        l[i] = std::pow(-1, i) * 0.2 * i;
        u[i] = 0.8 * i;
    }

    std::vector<double> x(n);
    //! ... copy/define the starting point
    for (int i = 0; i < n; i++) {
        x[i] = 1.;
    }

    //! ... we start the iteration by initializing task.
    strcpy(task, "START");

    double f;
    // gradient
    std::vector<double> g(n);
    //! ... the beginning of the optimization loop
    while (strncmp(task, "FG", 2) == 0 || strncmp(task, "NEW_X", 5) == 0 || strncmp(task, "START", 5) == 0) {
        // ... call the L-BFGS-B code.
        setulb_(&n, &m, x.data(), l.data(), u.data(), nbd.data(), &f, g.data(),
               &factr, &pgtol, wa.data(), iwa.data(), task, &iprint, csave, lsave,
                isave, dsave, (long int) 60, (long int) 60);

        std::cout << "task: " << task << "\n";
        
        // evaluate objective and gradient
        if (strncmp(task, "FG", 2) == 0) {
            std::cout << "COUCOU\n";
            f = problem.objective(x);
            g = problem.objective_dense_gradient(x);
        }
    }

    //! ... print output of this example
    std::cout << "Final L-BFGS-B Solution\n";
    std::cout << "lower bound   x-value      upper bound  gradient\n";
    double reduced_gradient = 0.;
    for (int i = 0; i < n - 1; i++) {
        std::cout << l[i] << "\t" << x[i] << "\t" << u[i] << "\t" << g[i] << "\n";
        reduced_gradient += std::abs(std::min(x[i] - l[i], u[i] - x[i]) * g[i]);
    }; // end for
    std::cout << "Reduced Gradient Norm = " << reduced_gradient << "\n";
    
    std::vector<double> constraint_multipliers;
    LocalSolution solution(x, g, constraint_multipliers);
    
    return solution;
}