#include "Problem.hpp"
#include <cmath>
#include <iostream>
#include "Utils.hpp"

Problem::Problem(std::string name, int number_variables, int number_constraints):
name(name), number_variables(number_variables), number_constraints(number_constraints),
// allocate all vectors
variable_name(number_variables), variable_discrete(number_variables), variables_bounds(number_variables), variable_status(number_variables),
constraint_name(number_constraints), constraint_variables(number_constraints), constraint_bounds(number_constraints), constraint_type(number_constraints), constraint_status(number_constraints)
{
}

Problem::~Problem() {
}

/* compute ||c|| */
double Problem::compute_constraint_residual(std::vector<double>& constraints, std::string norm_value) {
    std::vector<double> residuals(constraints.size());
    for (int j = 0; j < this->number_constraints; j++) {
        residuals[j] = std::max(std::max(0., this->constraint_bounds[j].lb - constraints[j]), constraints[j] - this->constraint_bounds[j].ub);
    }
    return norm(residuals, norm_value);
}

/* compute ||c_S|| for a given set S */
double Problem::compute_constraint_residual(std::vector<double>& constraints, std::set<int> constraint_set, std::string norm_value) {
    std::map<int, double> residuals;
    for (int j: constraint_set) {
        residuals[j] = std::max(std::max(0., this->constraint_bounds[j].lb - constraints[j]), constraints[j] - this->constraint_bounds[j].ub);
    }
    return norm(residuals, norm_value);
}

void Problem::determine_bounds_types(std::vector<Range>& bounds, std::vector<ConstraintType>& status) {
    for (unsigned int i = 0; i < bounds.size(); i++) {
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
    return;
}

void Problem::determine_constraints_() {
    int current_equality_constraint = 0;
    int current_inequality_constraint = 0;
    for (int j = 0; j < this->number_constraints; j++) {
        if (this->constraint_status[j] == EQUAL_BOUNDS) {
            this->equality_constraints[j] = current_equality_constraint;
            current_equality_constraint++;
        }
        else {
            this->inequality_constraints[j] = current_inequality_constraint;
            current_inequality_constraint++;
        }
    }
    return;
}
