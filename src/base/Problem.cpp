#include "Problem.hpp"
#include <cmath>
#include <iostream>
#include "Utils.hpp"

Problem::Problem(std::string name): name(name) {
}

Problem::~Problem() {
}

//! compute || c(feasible_constraints) ||_1 and || c(infeasible_constraints) ||_1 for index sets
// feasible_constraints, infeasible_constraints (overall length m)

double Problem::feasible_residual_norm(ConstraintPartition& constraint_partition, std::vector<double>& constraints, double chosen_norm) {
    std::vector<double> residuals(constraints.size());
    /* compute residuals for infeasible linear constraints */
    for (int j: constraint_partition.infeasible) {
        if (constraint_partition.constraint_feasibility[j] == INFEASIBLE_LOWER) {
            residuals[j] = std::max(0., constraints[j] - this->constraints_bounds[j].ub);
        }
        else {
            residuals[j] = std::max(0., this->constraints_bounds[j].lb - constraints[j]);
        }
    }
    /* compute residuals for feasible linear constraints */
    for (int j: constraint_partition.feasible) {
        residuals[j] = std::max(std::max(0., this->constraints_bounds[j].lb - constraints[j]), constraints[j] - this->constraints_bounds[j].ub);
    }
    return norm(residuals, chosen_norm);
}

double Problem::infeasible_residual_norm(ConstraintPartition& constraint_partition, std::vector<double>& constraints, double chosen_norm) {
    std::vector<double> residuals(constraints.size());
    /* compute residuals for infeasible constraints */
    for (int j: constraint_partition.infeasible) {
        residuals[j] = std::max(std::max(0., this->constraints_bounds[j].lb - constraints[j]), constraints[j] - this->constraints_bounds[j].ub);
    }
    return norm(residuals, chosen_norm);
}

/* compute ||c|| */
double Problem::infeasible_residual_norm(std::vector<double>& constraints, double chosen_norm) {
    std::vector<double> residuals(constraints.size());
    for (int j = 0; j < this->number_constraints; j++) {
        residuals[j] = std::max(std::max(0., this->constraints_bounds[j].lb - constraints[j]), constraints[j] - this->constraints_bounds[j].ub);
    }
    return norm(residuals, chosen_norm);
}

std::vector<ConstraintType> Problem::determine_bounds_types(std::vector<Range>& bounds) {
    std::vector<ConstraintType> status(bounds.size());

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
    return status;
}

void Problem::determine_constraints() {
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
