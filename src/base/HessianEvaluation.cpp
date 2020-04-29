#include <exception>
#include "HessianEvaluation.hpp"
#include "Utils.hpp"

HessianEvaluation::HessianEvaluation(int number_variables) : size(number_variables) {
}

HessianEvaluation::~HessianEvaluation() {
}

/* Exact Hessian */

ExactHessianEvaluation::ExactHessianEvaluation(int number_variables): HessianEvaluation(number_variables) {
}

void ExactHessianEvaluation::compute(Problem& problem, Iterate& iterate, double objective_multiplier, std::vector<double>& constraint_multipliers) {
    iterate.compute_hessian(problem, objective_multiplier, constraint_multipliers);
    return;
}

/* BFGS Hessian */

BFGSHessianEvaluation::BFGSHessianEvaluation(int number_variables): HessianEvaluation(number_variables), previous_x(number_variables) {
}

void BFGSHessianEvaluation::compute(Problem& problem, Iterate& iterate, double objective_multiplier, std::vector<double>& constraint_multipliers) {
    iterate.compute_hessian(problem, objective_multiplier, constraint_multipliers);
    return;
}