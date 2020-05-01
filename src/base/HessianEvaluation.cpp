#include <exception>
#include "HessianEvaluation.hpp"
#include "MA57Solver.hpp"
#include "Utils.hpp"

HessianEvaluation::HessianEvaluation(int number_variables) : size(number_variables), convexify(false), inertia(0.), inertia_last(0.) {
}

HessianEvaluation::~HessianEvaluation() {
}

CSCMatrix HessianEvaluation::modify_inertia(CSCMatrix& hessian) {
    MA57Solver solver;
    
    this->inertia = 0.;
    DEBUG << "Testing factorization with inertia term " << this->inertia << "\n";
    COOMatrix coo_hessian = hessian.to_COO();
    MA57Factorization factorization = solver.factorize(coo_hessian);

    bool good_inertia = false;
    if (!factorization.matrix_is_singular() && factorization.number_negative_eigenvalues() == 0) {
        DEBUG << "Factorization was a success\n";
        good_inertia = true;
    }
    else {
        // inertia term for Hessian
        if (this->inertia_last == 0.) {
            this->inertia = 1e-4;
        }
        else {
            this->inertia = std::max(1e-20, this->inertia_last / 3.);
        }
        hessian = hessian.add_identity_multiple(this->inertia);
        coo_hessian = hessian.to_COO();
    }

    while (!good_inertia) {
        DEBUG << "Testing factorization with inertia term " << this->inertia << "\n";
        factorization = solver.factorize(coo_hessian);

        if (!factorization.matrix_is_singular() && factorization.number_negative_eigenvalues() == 0) {
            good_inertia = true;
            DEBUG << "Factorization was a success\n";
            this->inertia_last = this->inertia;
        }
        else {
            if (this->inertia_last == 0.) {
                this->inertia *= 100.;
            }
            else {
                this->inertia *= 8.;
            }
            hessian = hessian.add_identity_multiple(this->inertia);
            coo_hessian = hessian.to_COO();
        }
    }
    DEBUG << "Final inertia: " << this->inertia << "\n";
    DEBUG << "Negative eigenvalues: " << factorization.number_negative_eigenvalues() << "\n";
            
    return hessian;
}

/* Exact Hessian */

ExactHessianEvaluation::ExactHessianEvaluation(int number_variables): HessianEvaluation(number_variables) {
}

void ExactHessianEvaluation::compute(Problem& problem, Iterate& iterate, double objective_multiplier, std::vector<double>& constraint_multipliers) {
    /* compute Hessian */
    iterate.compute_hessian(problem, objective_multiplier, constraint_multipliers);
    
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