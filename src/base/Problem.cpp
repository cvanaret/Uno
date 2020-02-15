#include "Problem.hpp"
#include <cmath>
#include <iostream>

Problem::Problem(std::string name): name(name) {
}

Problem::~Problem() {
}

//! compute || c(feasible_constraints) ||_1 and || c(infeasible_constraints) ||_1 for index sets
// feasible_constraints, infeasible_constraints (overall length m)
double Problem::feasible_residual_norm(ConstraintPartition& constraint_partition, std::vector<double>& constraints) {
	double feasible_residual = 0.;
	
	/* compute residuals for infeasible constraints */
	// TODO useful?
	for (unsigned int k = 0; k < constraint_partition.infeasible_set.size(); k++) {
		int j = constraint_partition.infeasible_set[k];
		if (constraint_partition.constraint_status[j] == INFEASIBLE_LOWER) {
			feasible_residual += std::max(0., constraints[j] - this->constraints_bounds[j].ub);
		}
		else {
			feasible_residual += std::max(0., this->constraints_bounds[j].lb - constraints[j]);
		}
	}

	/* compute residuals for feasible constraints */
	for (unsigned int k = 0; k < constraint_partition.feasible_set.size(); k++) {
		int j = constraint_partition.feasible_set[k];
		feasible_residual += std::max(0., this->constraints_bounds[j].lb - constraints[j]);
		feasible_residual += std::max(0., constraints[j] - this->constraints_bounds[j].ub);
	}
	return feasible_residual;
}

double Problem::infeasible_residual_norm(ConstraintPartition& constraint_partition, std::vector<double>& constraints) {
	double infeasible_residual = 0.;
	
	/* compute residuals for infeasible constraints */
	for (unsigned int k = 0; k < constraint_partition.infeasible_set.size(); k++) {
		int j = constraint_partition.infeasible_set[k];
		if (constraint_partition.constraint_status[j] == INFEASIBLE_LOWER) {
			infeasible_residual += std::max(0., this->constraints_bounds[j].lb - constraints[j]);
		}
		else {
			infeasible_residual += std::max(0., constraints[j] - this->constraints_bounds[j].ub);
		}
	}
	return infeasible_residual;
}

/* compute ||c||_1 */
double Problem::l1_norm(std::vector<double>& constraints) {
	double norm = 0.;

	for (int j = 0; j < this->number_constraints; j++) {
		double residual = std::max(constraints[j] - this->constraints_bounds[j].ub, this->constraints_bounds[j].lb - constraints[j]);
		norm += std::max(0., residual);
	}
	return norm;
}

std::vector<ConstraintType> Problem::determine_constraints_types(std::vector<Range>& variables_bounds) {
	std::vector<ConstraintType> status(variables_bounds.size());
	
	for (unsigned int i = 0; i < variables_bounds.size(); i++) {
		if (variables_bounds[i].lb == variables_bounds[i].ub) {
			status[i] = EQUAL_BOUNDS;
		}
		else if (-INFINITY < variables_bounds[i].lb && variables_bounds[i].ub < INFINITY) {
			status[i] = BOUNDED_BOTH_SIDES;
		}
		else if (-INFINITY < variables_bounds[i].lb) {
			status[i] = BOUNDED_LOWER;
		}
		else if (variables_bounds[i].ub < INFINITY) {
			status[i] = BOUNDED_UPPER;
		}
		else {
			status[i] = UNBOUNDED;
		}
	}
	return status;
}

std::map<int,int> Problem::determine_inequality_constraints() {
    std::map<int,int> inequality_constraints;
    int current_constraint = 0;
    for (int j = 0; j < this->number_constraints; j++) {
        if (this->constraint_status[j] != EQUAL_BOUNDS) {
            inequality_constraints[j] = current_constraint;
            current_constraint++;
        }
    }
    return inequality_constraints;
}
