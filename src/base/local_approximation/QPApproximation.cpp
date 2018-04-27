#include <cmath>
#include <map>
#include "QPApproximation.hpp"
#include "Constraint.hpp"
#include "Utils.hpp"
#include "Logger.hpp"

QPApproximation::QPApproximation(QPSolver& solver): LocalApproximation("QP"), solver(solver) {
}

void QPApproximation::allocate_solver(int number_variables, int number_constraints) {
	this->solver.allocate(number_variables, number_constraints);
}

LocalSolution QPApproximation::compute_optimality_step(Problem& problem, Iterate& current_point, double objective_multiplier, double radius) {
	/* generate the QP */
	QP qp = this->generate_optimality_qp(problem, current_point, objective_multiplier, radius);
	DEBUG << qp;
	
	/* generate the initial solution */
	std::vector<double> d0(qp.number_variables); // = {0.}
	
	/* solve the QP */
	LocalSolution solution = this->solver.solve(qp, d0);
	this->number_subproblems_solved++;
	return solution;
}

LocalSolution QPApproximation::compute_infeasibility_step(Problem& problem, Iterate& current_point, double radius,
		const std::vector<double>& d, ConstraintPartition& constraint_partition, std::vector<double>& multipliers) {
	/* generate the QP */
	QP qp = this->generate_infeasibility_qp(problem, current_point, radius, constraint_partition, multipliers);
	DEBUG << qp;
	
	/* generate the initial solution */
	std::vector<double> d0 = d;
	
	/* solve the QP */
	LocalSolution solution = this->solver.solve(qp, d0);
	this->number_subproblems_solved++;
	return solution;
}

/* additional variables */
LocalSolution QPApproximation::compute_l1_penalty_step(Problem& problem, Iterate& current_point, double radius, double penalty_parameter, PenaltyConstraints penalty_constraints) {
	/* generate the QP */
	QP qp = this->generate_l1_penalty_qp(problem, current_point, radius, penalty_parameter, penalty_constraints);
	DEBUG << qp;
	
	/* generate the initial solution */
	std::vector<double> d0(qp.number_variables); // = {0.}

	/* solve the QP */
	LocalSolution solution = this->solver.solve(qp, d0);
	this->number_subproblems_solved++;
	return solution;
}

QP QPApproximation::generate_qp(Problem& problem, Iterate& current_point, double radius) {
	/* initialize the QP */
	QP qp(problem.number_variables, problem.number_constraints, current_point.hessian);
	
	/* bound constraints intersected with trust region  */
	for (int i = 0; i < problem.number_variables; i++) {
		qp.variable_lb[i] = std::max(-radius, problem.variable_lb[i] - current_point.x[i]);
		qp.variable_ub[i] = std::min(radius, problem.variable_ub[i] - current_point.x[i]);
	}
	
	/* compute the constraints */
	this->set_constraints(problem, qp, current_point);
	
	return qp;
}

QP QPApproximation::generate_optimality_qp(Problem& problem, Iterate& current_point, double objective_multiplier, double radius) {
	DEBUG << "Creating the optimality problem\n";
	
	/* compute the Lagrangian Hessian */
	/* keep only the multipliers of the general constraints */
	std::vector<double> constraint_multipliers(problem.number_constraints);
	for (int j = 0; j < problem.number_constraints; j++) {
		/* take the opposite of the multiplier. The reason is that in AMPL, the
		 * Lagrangian is f + lambda.g, while Argonot uses f - lambda.g */
		constraint_multipliers[j] = -current_point.multipliers[problem.number_variables + j];
	}
	current_point.compute_hessian(problem, objective_multiplier, constraint_multipliers);
	
	/* initialize the QP */
	QP qp = this->generate_qp(problem, current_point, radius);
	
	/* bounds of the linearized constraints */
	for (int j = 0; j < qp.number_constraints; j++) {
		qp.constraint_lb[j] = problem.constraint_lb[j] - current_point.constraints[j];
		qp.constraint_ub[j] = problem.constraint_ub[j] - current_point.constraints[j];
	}
	
	/* compute the objective */
	this->set_optimality_objective(problem, qp, current_point);
	
	return qp;
}

QP QPApproximation::generate_infeasibility_qp(Problem& problem, Iterate& current_point, double radius, ConstraintPartition& constraint_partition, std::vector<double>& multipliers) {
	int number_infeasible = constraint_partition.infeasible_set.size();
	DEBUG << "Creating the restoration problem with infeasible " << number_infeasible << " constraints\n";

	/* compute the Lagrangian Hessian */
	double objective_multiplier = 0.;
	/* keep only the multipliers of the general constraints */
	std::vector<double> constraint_multipliers(problem.number_constraints);
	for (int j = 0; j < problem.number_constraints; j++) {
		/* take the opposite of the multiplier. The reason is that in AMPL, the
		 * Lagrangian is f + lambda.g, while Argonot uses f - lambda.g */
		if (constraint_partition.status[j] == INFEASIBLE_LOWER) {
			constraint_multipliers[j] = -1.;
		}
		else if (constraint_partition.status[j] == INFEASIBLE_UPPER) {
			constraint_multipliers[j] = 1.;
		}
		else {
			constraint_multipliers[j] = -current_point.multipliers[problem.number_variables + j];
		}
	}
	current_point.compute_hessian(problem, objective_multiplier, constraint_multipliers);
	
	/* initialize the QP */
	QP qp = this->generate_qp(problem, current_point, radius);
	
	/* bounds of the linearized constraints */
	for (int j = 0; j < qp.number_constraints; j++) {
		if (constraint_partition.status[j] == INFEASIBLE_LOWER) {
			qp.constraint_lb[j] = -INFINITY;
			qp.constraint_ub[j] = problem.constraint_lb[j] - current_point.constraints[j];
		}
		else if (constraint_partition.status[j] == INFEASIBLE_UPPER) {
			qp.constraint_lb[j] = problem.constraint_ub[j] - current_point.constraints[j];
			qp.constraint_ub[j] = INFINITY;
		}
		else { // FEASIBLE
			qp.constraint_lb[j] = problem.constraint_lb[j] - current_point.constraints[j];
			qp.constraint_ub[j] = problem.constraint_ub[j] - current_point.constraints[j];
		}
	}
	
	/* compute the objective */
	this->set_infeasibility_objective(problem, qp, current_point, constraint_partition);
	return qp;
}

QP QPApproximation::generate_l1_penalty_qp(Problem& problem, Iterate& current_point, double radius, double penalty_parameter, PenaltyConstraints penalty_constraints) {
	int number_variables = problem.number_variables + penalty_constraints.number_additional_variables;
	int number_constraints = penalty_constraints.number_constraints;
	
	/* compute the Lagrangian Hessian from scratch */
	current_point.is_hessian_computed = false;
	double objective_multiplier = penalty_parameter;
	/* keep only the multipliers of the general constraints */
	std::vector<double> constraint_multipliers(problem.number_constraints);
	for (int j = 0; j < problem.number_constraints; j++) {
		/* take the opposite of the multiplier. The reason is that in AMPL, the
		 * Lagrangian is f + lambda.g, while Argonot uses f - lambda.g */
		constraint_multipliers[j] = -current_point.multipliers[problem.number_variables + j];
	}
	current_point.compute_hessian(problem, objective_multiplier, constraint_multipliers);
	
	/* initialize the QP */
	QP qp(number_variables, number_constraints, current_point.hessian);
	
	/* bounds of original variables intersected with trust region  */
	for (int i = 0; i < problem.number_variables; i++) {
		qp.variable_lb[i] = std::max(-radius, problem.variable_lb[i] - current_point.x[i]);
		qp.variable_ub[i] = std::min(radius, problem.variable_ub[i] - current_point.x[i]);
	}
	/* bounds of additional variables */
	for (int k = 0; k < penalty_constraints.number_additional_variables; k++) {
		qp.variable_lb[problem.number_variables + k] = 0.;
		qp.variable_ub[problem.number_variables + k] = INFINITY;
	}
	
	/* apply the nonzero penalty parameter on the initial objective */
	if (penalty_parameter != 0.) {
		if (!current_point.is_objective_gradient_computed) {
			std::map<int,double> objective_gradient = problem.objective_sparse_gradient(current_point.x);
			current_point.set_objective_gradient(objective_gradient);
		}
		qp.objective = current_point.objective_gradient;
		for (std::map<int,double>::iterator it = qp.objective.begin(); it != qp.objective.end(); it++) {
			int index = it->first;
			qp.objective[index] *= penalty_parameter;
		}
	}
	/* add additional variables to the objective */
	for (int k = 0; k < penalty_constraints.number_additional_variables; k++) {
		qp.objective[problem.number_variables + k] = 1.;
	}
	
	/* compute the original constraint gradients */
	current_point.compute_constraint_jacobian(problem);

	/* add the constraints */
	int current_additional_variable = problem.number_variables;
	int current_constraint = 0;
	for (int j = 0; j < problem.number_constraints; j++) {
		if (penalty_constraints.status[j] == EQUALITY) {
			/* a single constraint with both additional variables */
			std::map<int,double> gradient(current_point.constraint_jacobian[j]);
			gradient[current_additional_variable] = -1.;
			gradient[current_additional_variable+1] = 1.;
			qp.constraints[current_constraint] = gradient;
			/* identical bounds */
			qp.constraint_lb[current_constraint] = problem.constraint_lb[j] - current_point.constraints[j];
			qp.constraint_ub[current_constraint] = problem.constraint_ub[j] - current_point.constraints[j];
			current_additional_variable += 2;
			current_constraint++;
		}
		if (penalty_constraints.status[j] == BOUNDED_BOTH_SIDES || penalty_constraints.status[j] == BOUNDED_LOWER) {
			/* a single constraint with one additional variable */
			std::map<int,double> gradient(current_point.constraint_jacobian[j]);
			gradient[current_additional_variable] = 1.;
			qp.constraints[current_constraint] = gradient;
			/* bounds */
			qp.constraint_lb[current_constraint] = problem.constraint_lb[j] - current_point.constraints[j];
			qp.constraint_ub[current_constraint] = INFINITY;
			current_additional_variable++;
			current_constraint++;
		}
		if (penalty_constraints.status[j] == BOUNDED_BOTH_SIDES || penalty_constraints.status[j] == BOUNDED_UPPER) {
			/* a single constraint with one additional variable */
			std::map<int,double> gradient(current_point.constraint_jacobian[j]);
			gradient[current_additional_variable] = -1.;
			qp.constraints[current_constraint] = gradient;
			/* bounds */
			qp.constraint_lb[current_constraint] = -INFINITY;
			qp.constraint_ub[current_constraint] = problem.constraint_ub[j] - current_point.constraints[j];
			current_additional_variable++;
			current_constraint++;
		}
	}
	return qp;
}

void QPApproximation::set_constraints(Problem& problem, QP& qp, Iterate& current_point) {
	/* compute the constraint Jacobian */
	current_point.compute_constraint_jacobian(problem);
	qp.constraints = current_point.constraint_jacobian;
	return;
}

void QPApproximation::set_optimality_objective(Problem& problem, QP& qp, Iterate& current_point) {
	/* compute the objective Jacobian */
	if (!current_point.is_objective_gradient_computed) {
		std::map<int,double> objective_gradient = problem.objective_sparse_gradient(current_point.x);
		current_point.set_objective_gradient(objective_gradient);
	}	
	qp.objective = current_point.objective_gradient;
	return;
}

void QPApproximation::set_infeasibility_objective(Problem& problem, QP& qp, Iterate& current_point, ConstraintPartition& constraint_partition) {
	/* objective function: add the gradients of infeasible constraints */
	std::map<int,double> objective_gradient;
	
	for (unsigned int k = 0; k < constraint_partition.infeasible_set.size(); k++) {
		int j = constraint_partition.infeasible_set[k];
		/* combine into objective_gradient */
		
		for (std::map<int,double>::iterator it = current_point.constraint_jacobian[j].begin(); it != current_point.constraint_jacobian[j].end(); it++) {
			int variable_index = it->first;
			double derivative = it->second;
			
			if (constraint_partition.status[j] == INFEASIBLE_LOWER) {
				objective_gradient[variable_index] -= derivative;
			}
			else {
				objective_gradient[variable_index] += derivative;
			}
		}
	}
	current_point.set_objective_gradient(objective_gradient);
	qp.objective = current_point.objective_gradient;
	return;
}
