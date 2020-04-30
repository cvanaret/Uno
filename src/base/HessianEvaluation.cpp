#include <exception>
#include "HessianEvaluation.hpp"
#include "Utils.hpp"

HessianEvaluation::HessianEvaluation(int number_variables) : size(number_variables), convexify(false) {
}

HessianEvaluation::~HessianEvaluation() {
}

/* Exact Hessian */

ExactHessianEvaluation::ExactHessianEvaluation(int number_variables): HessianEvaluation(number_variables) {
}

void ExactHessianEvaluation::compute(Problem& problem, Iterate& iterate, double objective_multiplier, std::vector<double>& constraint_multipliers) {
    //if (this->convexify) {
    //    throw std::runtime_error("ExactHessianEvaluation::compute should convexify");
    //}
    //else {
        iterate.compute_hessian(problem, objective_multiplier, constraint_multipliers);
    //}
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