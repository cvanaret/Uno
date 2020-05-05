#include <exception>
#include "HessianEvaluation.hpp"
#include "MA57Solver.hpp"
#include "Utils.hpp"

HessianEvaluation::HessianEvaluation(int number_variables) : size(number_variables), convexify(false) {
}

HessianEvaluation::~HessianEvaluation() {
}

CSCMatrix HessianEvaluation::modify_inertia(CSCMatrix& hessian) {
    MA57Solver solver;
    double beta = 1e-4;
    
    // Nocedal and Wright, p51
    double smallest_diagonal_entry = hessian.smallest_diagonal_entry();
    DEBUG << "The minimal diagonal entry of the Hessian is " << hessian.smallest_diagonal_entry() << "\n";
    
    double inertia = 0.;
    double previous_inertia = 0.;
    if (smallest_diagonal_entry <= 0.) {
        inertia = beta - smallest_diagonal_entry;
    }
    
    COOMatrix coo_hessian = hessian.to_COO();
    MA57Factorization factorization = solver.factorize(coo_hessian);
    
    bool good_inertia = false;
    while (!good_inertia) {
        DEBUG << "Testing factorization with inertia term " << inertia << "\n";
        
        if (!factorization.matrix_is_singular() && factorization.number_negative_eigenvalues() == 0) {
            good_inertia = true;
            DEBUG << "Factorization was a success\n";
        }
        else {
            if (inertia == 0.) {
                inertia = beta;
                previous_inertia = 0.;
            }
            else {
                previous_inertia = inertia;
                inertia *= 10.;
            }
            hessian = hessian.add_identity_multiple(inertia - previous_inertia);
            coo_hessian = hessian.to_COO();
            factorization = solver.factorize(coo_hessian);
        }
    }
    DEBUG << "Final inertia: " << inertia << "\n";
    DEBUG << "Negative eigenvalues: " << factorization.number_negative_eigenvalues() << "\n";
    return hessian;
}

/* Exact Hessian */

ExactHessianEvaluation::ExactHessianEvaluation(int number_variables): HessianEvaluation(number_variables) {
}

void ExactHessianEvaluation::compute(Problem& problem, Iterate& iterate, double objective_multiplier, std::vector<double>& constraint_multipliers) {
    /* compute Hessian */
    iterate.compute_hessian(problem, objective_multiplier, constraint_multipliers);
    
    DEBUG << "hessian before convexification: " << iterate.hessian;
    
    if (this->convexify) {
        /* modify the inerta to make*/
        iterate.hessian = this->modify_inertia(iterate.hessian);
    }
    return;
}

/* BFGS Hessian */

BFGSHessianEvaluation::BFGSHessianEvaluation(int number_variables): HessianEvaluation(number_variables), previous_x(number_variables) {
}

void BFGSHessianEvaluation::compute(Problem& problem, Iterate& iterate, double objective_multiplier, std::vector<double>& constraint_multipliers) {
    // the BFGS Hessian is already positive definite, do not convexify
    iterate.compute_hessian(problem, objective_multiplier, constraint_multipliers);
    return;
}