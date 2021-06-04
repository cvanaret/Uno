#include <cmath>
#include <iostream>
#include <cassert>
#include "Problem.hpp"
#include "Vector.hpp"

std::map<FunctionType, std::string> Problem::type_to_string = {
    {LINEAR, "linear"},
    {QUADRATIC, "quadratic"},
    {NONLINEAR, "nonlinear"}
};

/* Abstract Problem class */

Problem::Problem(std::string& name, int number_variables, int number_constraints, FunctionType type):
name(name), number_variables(number_variables), number_constraints(number_constraints), type(type),
objective_sign(1.), objective_type(NONLINEAR),
// allocate all vectors
variable_name(number_variables),
//variable_discrete(number_variables),
variables_bounds(number_variables), variable_status(number_variables),
constraint_name(number_constraints),
//constraint_variables(number_constraints),
constraint_bounds(number_constraints), constraint_type(number_constraints), constraint_status(number_constraints),
hessian_maximum_number_nonzeros(0) {
}

/* compute ||c|| */
double Problem::compute_constraint_residual(const std::vector<double>& constraints, Norm residual_norm) const {
    std::vector<double> residuals(constraints.size());
    for (size_t j = 0; j < this->number_constraints; j++) {
        residuals[j] = std::max(std::max(0., this->constraint_bounds[j].lb - constraints[j]), constraints[j] - this->constraint_bounds[j].ub);
    }
    return norm(residuals, residual_norm);
}

/* compute ||c_S|| for a given set S */
double Problem::compute_constraint_residual(const std::vector<double>& constraints, const std::set<int>& constraint_set, Norm residual_norm) const {
    SparseGradient residuals;
    for (int j: constraint_set) {
        residuals[j] = std::max(std::max(0., this->constraint_bounds[j].lb - constraints[j]), constraints[j] - this->constraint_bounds[j].ub);
    }
    return norm(residuals, residual_norm);
}

void Problem::determine_bounds_types(std::vector<Range>& bounds, std::vector<ConstraintType>& status) {
   assert(bounds.size() == status.size());

    for (size_t i = 0; i < bounds.size(); i++) {
        if (bounds[i].lb == bounds[i].ub) {
            status[i] = EQUAL_BOUNDS;
        }
        else if (-INFINITY < bounds[i].lb && bounds[i].ub < INFINITY) {
            status[i] = BOUNDED_BOTH_SIDES;
        }
        else if (-INFINITY < bounds[i].lb) {
            status[i] = BOUNDED_LOWER;
        }
        else if (bounds[i].ub < INFINITY) {
            status[i] = BOUNDED_UPPER;
        }
        else {
            status[i] = UNBOUNDED;
        }
    }
}

void Problem::determine_constraints_() {
    int current_equality_constraint = 0;
    int current_inequality_constraint = 0;
    for (size_t j = 0; j < this->number_constraints; j++) {
        if (this->constraint_status[j] == EQUAL_BOUNDS) {
            this->equality_constraints[j] = current_equality_constraint;
            current_equality_constraint++;
        }
        else {
            this->inequality_constraints[j] = current_inequality_constraint;
            current_inequality_constraint++;
        }
    }
}

/* native C++ problem */

//CppProblem::CppProblem(std::string name, int number_variables, int number_constraints, double (*objective)(std::vector<double> x), std::vector<double> (*objective_gradient)(std::vector<double> x)):
//Problem(name, number_variables, number_constraints),
//objective_(objective),
//objective_gradient_(objective_gradient) {
//}
//
//double CppProblem::objective(std::vector<double>& x) {
//    return this->objective_(x);
//}
//
//std::vector<double> CppProblem::objective_dense_gradient(std::vector<double>& x) {
//    return this->objective_gradient_(x);
//}
//
//SparseGradient CppProblem::objective_sparse_gradient(std::vector<double>& x) {
//    std::vector<double> dense_gradient = this->objective_gradient_(x);
//    SparseGradient sparse_gradient;
//    for (size_t i = 0; i < dense_gradient.size(); i++) {
//        if (dense_gradient[i] != 0.) {
//            sparse_gradient[i] = dense_gradient[i];
//        }
//    }
//    return sparse_gradient;
//}
//
//double CppProblem::evaluate_constraint(int j, std::vector<double>& x) {
//    return this->constraints_[j](x);
//}
//
//std::vector<double> CppProblem::evaluate_constraints(std::vector<double>& x) {
//    std::vector<double> constraints(this->number_constraints);
//    for (int j = 0; j < this->number_constraints; j++) {
//        constraints[j] = this->evaluate_constraint(j, x);
//    }
//    return constraints;
//}
