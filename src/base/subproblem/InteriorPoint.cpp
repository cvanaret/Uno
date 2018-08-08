#include <cmath>
#include "InteriorPoint.hpp"

InteriorPoint::InteriorPoint(): Subproblem("IPM"), mu(0.), tau(0.995), default_multiplier(0.5), default_slack(1.) {
}

void InteriorPoint::initialize(Problem& problem, Iterate& current_iterate, int number_variables, int number_constraints, double radius) {
	this->variable_status.resize(problem.number_variables);
	
	/* determine the type of bound constraints (equality, bounded below or above) */
	for (int i = 0; i < problem.number_variables; i++) {
		if (problem.variable_lb[i] == problem.variable_ub[i]) {
			this->variable_status[i] = EQUAL_BOUNDS;
		}
		else if (radius < INFINITY || (-INFINITY < problem.variable_lb[i] && problem.variable_ub[i] < INFINITY)) {
			this->variable_status[i] = BOUNDED_BOTH_SIDES;
		}
		else if (-INFINITY < problem.variable_lb[i]) {
			this->variable_status[i] = BOUNDED_LOWER;
		}
		else if (problem.variable_ub[i] < INFINITY) {
			this->variable_status[i] = BOUNDED_UPPER;
		}
		else {
			this->variable_status[i] = UNBOUNDED;
		}
	}
	
	std::vector<double> variable_lb(problem.number_variables);
	std::vector<double> variable_ub(problem.number_variables);
	for (int i = 0; i < problem.number_variables; i++) {
		variable_lb[i] = std::max(current_iterate.x[i] - radius, problem.variable_lb[i]);
		variable_ub[i] = std::min(current_iterate.x[i] + radius, problem.variable_ub[i]);
	}
	
	for (int i = 0; i < problem.number_variables; i++) {
		std::cout << "x" << i << " in [" << variable_lb[i] << ", " << variable_ub[i] << "]\n";
	}
	std::cout << "\n";
	
	/* determine the slacks and multipliers */
	for (int i = 0; i < problem.number_variables; i++) {
		if (this->variable_status[i] == EQUAL_BOUNDS) {
			std::cout << "x" << i << ": EQUAL_BOUNDS\n";
			throw std::out_of_range("IPM: variable has equal bounds, behavior not specified");
		}
		if (this->variable_status[i] == BOUNDED_LOWER || this->variable_status[i] == BOUNDED_BOTH_SIDES) {
			std::cout << "x" << i << ": BOUNDED_LOWER\n";
			this->multipliers.push_back(this->default_multiplier); // positive multiplier
			this->slacks.push_back(current_iterate.x[i] - variable_lb[i]); // exact slack
		}
		if (this->variable_status[i] == BOUNDED_UPPER || this->variable_status[i] == BOUNDED_BOTH_SIDES) {
			std::cout << "x" << i << ": BOUNDED_UPPER\n";
			this->multipliers.push_back(-this->default_multiplier); // negative multiplier
			this->slacks.push_back(variable_ub[i] - current_iterate.x[i]); // exact slack
		}
	}
	for (int j = 0; j < problem.number_constraints; j++) {
		if (problem.constraint_status[j] == EQUAL_BOUNDS) {
			std::cout << "c" << j << ": EQUAL_BOUNDS\n";
			this->multipliers.push_back(this->default_multiplier); // unconstrained multiplier
			// no slack
		}
		/* slacks for inequality constraints */
		if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			std::cout << "c" << j << ": BOUNDED_LOWER\n";
			this->multipliers.push_back(this->default_multiplier); // positive multiplier
			this->slacks.push_back(this->default_slack);
		}
		if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			std::cout << "c" << j << ": BOUNDED_UPPER\n";
			this->multipliers.push_back(-this->default_multiplier); // negative multiplier
			this->slacks.push_back(this->default_slack);
		}
	}
	
	std::cout << this->slacks.size() << " slack variables\n";
	std::cout << this->multipliers.size() << " constraints/multipliers\n";
	
	/* compute a step with mu = 0 to initialize slacks, multipliers and mu */
	this->compute_optimality_step(problem, current_iterate, 1., radius);
}
		
LocalSolution InteriorPoint::compute_optimality_step(Problem& problem, Iterate& current_iterate, double objective_multiplier, double radius) {
	/* reduced primal-dual approach */
	
	/* build trust region */
	std::vector<double> variable_lb(problem.number_variables);
	std::vector<double> variable_ub(problem.number_variables);
	for (int i = 0; i < problem.number_variables; i++) {
		variable_lb[i] = std::max(current_iterate.x[i] - radius, problem.variable_lb[i]);
		variable_ub[i] = std::min(current_iterate.x[i] + radius, problem.variable_ub[i]);
	}
	for (int i = 0; i < problem.number_variables; i++) {
		std::cout << "x" << i << " in [" << variable_lb[i] << ", " << variable_ub[i] << "]\n";
	}
	std::cout << "\n";

	std::cout << "x is: "; print_vector(std::cout, current_iterate.x);
	std::cout << "s is: "; print_vector(std::cout, slacks);
	std::cout << "λ is: "; print_vector(std::cout, multipliers);
	std::cout << "mu is " << this->mu << "\n\n";
	
	/* compute first-order information */
	current_iterate.compute_objective_gradient(problem);
	current_iterate.compute_constraints_jacobian(problem);
	
	/*******************************/
	/* sparse symmetric KKT matrix */
	/*******************************/
	
	/* evaluate the Lagrangian Hessian */
	// recreate the multipliers
	std::vector<double> original_multipliers(problem.number_constraints);
	int current_multiplier = 0;
	for (int i = 0; i < problem.number_variables; i++) {
		if (this->variable_status[i] == BOUNDED_LOWER || this->variable_status[i] == BOUNDED_BOTH_SIDES) {
			current_multiplier++;
		}
		if (this->variable_status[i] == BOUNDED_UPPER || this->variable_status[i] == BOUNDED_BOTH_SIDES) {
			current_multiplier++;
		}
	}
	for (int j = 0; j < problem.number_constraints; j++) {
		original_multipliers[j] = multipliers[current_multiplier];
		current_multiplier++;
		// if both bounds, add the two multipliers
		if (problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			original_multipliers[j] += multipliers[current_multiplier];
			current_multiplier++;
		}
	}
	// create the KKT matrix
	COOMatrix kkt_matrix = this->generate_kkt_matrix(problem, current_iterate, objective_multiplier, original_multipliers);
	
	/*******************/
	/* right-hand side */
	/*******************/
	
	std::vector<double> rhs = this->generate_rhs(problem, current_iterate, objective_multiplier, original_multipliers, variable_lb, variable_ub);
	
	std::cout << "KKT matrix:\n";
	for (unsigned int k = 0; k < kkt_matrix.matrix.size(); k++) {
		std::cout << "m(" << kkt_matrix.row_indices[k] << ", " << kkt_matrix.column_indices[k] << ") = " << kkt_matrix.matrix[k] << "\n";
	}
	print_vector(std::cout, rhs);
	
	/************/
	/* solution */
	/************/
	
	/* compute the solution (Δx, -Δλ) */
	LocalSolution solution = this->solver.solve(kkt_matrix, rhs);
	
	/* retrieve +Δλ (Nocedal p590) */
	for (unsigned int j = problem.number_variables; j < problem.number_variables + multipliers.size(); j++) {
		solution.x[j] = -solution.x[j];
	}
	
	std::cout << "MA57 solution:\n";
	for (int k = 0; k < problem.number_variables; k++) {
		std::cout << "Δx" << k << " = " << solution.x[k] << "\n";
	}
	for (unsigned int j = 0; j < multipliers.size(); j++) {
		std::cout << "Δλ" << j << " = " << solution.x[problem.number_variables + j] << "\n";
	}
	
	/* compute slack displacements Δs */
	std::vector<double> delta_s = this->compute_slack_displacements(problem, mu, solution.x);

	for (unsigned int j = 0; j < slacks.size(); j++) {
		std::cout << "Δs" << j << " = " << delta_s[j] << "\n";
	}
	std::cout << "\n";
	
	/* fraction to boundary rule for slacks */
	double alpha_s = 1.;
	for (unsigned int j = 0; j < slacks.size(); j++) {
		double trial_alpha_sj = -this->tau*slacks[j]/delta_s[j];
		if (0 < trial_alpha_sj && trial_alpha_sj <= 1.) {
			alpha_s = std::min(alpha_s, trial_alpha_sj);
		}
	}
	std::cout << "alpha_s = " << alpha_s << "\n";
	
	/* fraction to boundary rule for multipliers */
	double alpha_mult = 1.;
	current_multiplier = 0;
	for (int i = 0; i < problem.number_variables; i++) {
		if (variable_status[i] == BOUNDED_LOWER || variable_status[i] == BOUNDED_BOTH_SIDES) {
			double trial_alpha_lj = -this->tau*multipliers[current_multiplier]/solution.x[problem.number_variables + current_multiplier];
			if (0 < trial_alpha_lj && trial_alpha_lj <= 1.) {
				alpha_mult = std::min(alpha_mult, trial_alpha_lj);
			}
			current_multiplier++;
		}
		if (variable_status[i] == BOUNDED_UPPER || variable_status[i] == BOUNDED_BOTH_SIDES) {
			double trial_alpha_lj = -this->tau*multipliers[current_multiplier]/solution.x[problem.number_variables + current_multiplier];
			if (0 < trial_alpha_lj && trial_alpha_lj <= 1.) {
				alpha_mult = std::min(alpha_mult, trial_alpha_lj);
			}
			current_multiplier++;
		}
	}
	for (int j = 0; j < problem.number_constraints; j++) {
		if (problem.constraint_status[j] == EQUAL_BOUNDS) {
			// multiplier is unconstrained
			current_multiplier++;
		}
		else {
			int number_copies = (variable_status[j] == BOUNDED_BOTH_SIDES) ? 2 : 1;
			for (int copy = 1; copy <= number_copies; copy++) {
				double trial_alpha_lj = -this->tau*multipliers[current_multiplier]/solution.x[problem.number_variables + current_multiplier];
				if (0 < trial_alpha_lj && trial_alpha_lj <= 1.) {
					alpha_mult = std::min(alpha_mult, trial_alpha_lj);
				}
				current_multiplier++;
			}
		}
	}
	std::cout << "alpha_λ = " << alpha_mult << "\n\n";

	std::cout << "New x: (";
	std::vector<double> trial_x(problem.number_variables);
	for (int i = 0; i < problem.number_variables; i++) {
		trial_x[i] = current_iterate.x[i] + alpha_s*solution.x[i];
		std::cout << trial_x[i] << ", ";
	}
	std::cout << ")\n";
	std::cout << "New s: (";
	for (unsigned int j = 0; j < slacks.size(); j++) {
		this->slacks[j] += alpha_s*delta_s[j];
		std::cout << this->slacks[j] << ", ";
	}
	std::cout << ")\n";
	std::cout << "New λ: (";
	for (unsigned int j = 0; j < multipliers.size(); j++) {
		this->multipliers[j] += alpha_mult*solution.x[problem.number_variables + j];
		std::cout << this->multipliers[j] << ", ";
	}
	std::cout << ")\n";
	
	/*******************************/
	/* update of barrier parameter */
	/*******************************/
	
	double sum = 0.;
	current_multiplier = 0;
	int current_slack = 0;
	
	for (int i = 0; i < problem.number_variables; i++) {
		if (variable_status[i] == BOUNDED_LOWER || variable_status[i] == BOUNDED_BOTH_SIDES) {
			sum += multipliers[current_multiplier]*slacks[current_slack];
			current_multiplier++;
			current_slack++;
		}
		if (variable_status[i] == BOUNDED_UPPER || variable_status[i] == BOUNDED_BOTH_SIDES) {
			sum -= multipliers[current_multiplier]*slacks[current_slack];
			current_multiplier++;
			current_slack++;
		}
	}
	for (int j = 0; j < problem.number_constraints; j++) {
		if (problem.constraint_status[j] == EQUAL_BOUNDS) {
			current_multiplier++;
		}
		/* slacks for inequality constraints */
		if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			sum += multipliers[current_multiplier]*slacks[current_slack];
			current_multiplier++;
			current_slack++;
		}
		if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			sum -= multipliers[current_multiplier]*slacks[current_slack];
			current_multiplier++;
			current_slack++;
		}
	}
	
	sum /= slacks.size();
	std::cout << "NEW mu VALUE: " << sum << "\n";
	this->mu = sum;
	
	LocalSolution sol(trial_x, 0, 0);
	return sol;
}

LocalSolution InteriorPoint::compute_infeasibility_step(Problem& problem, Iterate& current_iterate, double radius, const std::vector<double>& d, ConstraintPartition& constraint_partition, std::vector<double>& multipliers) {
	std::vector<double> x;
	LocalSolution solution(x, 0, 0);
	return solution;
}

LocalSolution InteriorPoint::compute_l1_penalty_step(Problem& problem, Iterate& current_iterate, double radius, double penalty_parameter, PenaltyDimensions penalty_dimensions) {
	std::vector<double> x;
	LocalSolution solution(x, 0, 0);
	return solution;
}

COOMatrix InteriorPoint::generate_kkt_matrix(Problem& problem, Iterate& current_iterate, double objective_multiplier, std::vector<double> original_multipliers) {
	current_iterate.compute_hessian(problem, objective_multiplier, original_multipliers);
	COOMatrix kkt_matrix = current_iterate.hessian.to_COO();
	kkt_matrix.size += this->multipliers.size();
	int current_column = problem.number_variables;
	
	/* Jacobian of bound constraints */
	for (int i = 0; i < problem.number_variables; i++) {
		if (this->variable_status[i] != UNBOUNDED) {
			kkt_matrix.add_term(1., i, current_column);
			current_column++;
		}
		if (this->variable_status[i] == BOUNDED_BOTH_SIDES) {
			kkt_matrix.add_term(1., i, current_column);
			current_column++;
		}
	}
	/* Jacobian of general constraints */
	for (int j = 0; j < problem.number_constraints; j++) {
		int number_copies = (problem.constraint_status[j] == BOUNDED_BOTH_SIDES) ? 2 : 1;
		for (int copy = 1; copy <= number_copies; copy++) {
			for (std::map<int,double>::iterator it = current_iterate.constraints_jacobian[j].begin(); it != current_iterate.constraints_jacobian[j].end(); it++) {
				int variable_index = it->first;
				double derivative = it->second;
				kkt_matrix.add_term(derivative, variable_index, current_column);
			}
			current_column++;
		}
	}
	
	/* variable diagonal terms +/-L^-1 S */
	int current_slack = 0;
	int current_multiplier = 0;
	for (int i = 0; i < problem.number_variables; i++) {
		if (this->variable_status[i] == BOUNDED_LOWER || this->variable_status[i] == BOUNDED_BOTH_SIDES) {
			current_column = problem.number_variables + current_multiplier;
			kkt_matrix.add_term(-this->slacks[current_slack]/this->multipliers[current_multiplier], current_column, current_column);
			current_slack++;
			current_multiplier++;
		}
		if (this->variable_status[i] == BOUNDED_UPPER || this->variable_status[i] == BOUNDED_BOTH_SIDES) {
			current_column = problem.number_variables + current_multiplier;
			kkt_matrix.add_term(this->slacks[current_slack]/this->multipliers[current_multiplier], current_column, current_column);
			current_slack++;
			current_multiplier++;
		}
	}
	/* constraint diagonal terms +/-L^-1 S */
	for (int j = 0; j < problem.number_constraints; j++) {
		if (problem.constraint_status[j] == EQUAL_BOUNDS) {
			// no slack
			current_multiplier++;
		}
		else {
			if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
				current_column = problem.number_variables + current_multiplier;
				kkt_matrix.add_term(-this->slacks[current_slack]/this->multipliers[current_multiplier], current_column, current_column);
				current_slack++;
				current_multiplier++;
			}
			if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
				current_column = problem.number_variables + current_multiplier;
				kkt_matrix.add_term(this->slacks[current_slack]/this->multipliers[current_multiplier], current_column, current_column);
				current_slack++;
				current_multiplier++;
			}
		}
	}
	return kkt_matrix;
}

std::vector<double> InteriorPoint::generate_rhs(Problem& problem, Iterate& current_iterate, double objective_multiplier, std::vector<double>& original_multipliers,
		std::vector<double> variable_lb, std::vector<double> variable_ub) {
	std::vector<double> rhs(problem.number_variables + this->multipliers.size());
	
	/* objective gradient */
	for (std::map<int,double>::iterator it = current_iterate.objective_gradient.begin(); it != current_iterate.objective_gradient.end(); it++) {
		int variable_index = it->first;
		double derivative = it->second;
		rhs[variable_index] -= objective_multiplier*derivative;
	}
	/* constraint gradients */
	for (int j = 0; j < problem.number_constraints; j++) {
		double multiplier_j = original_multipliers[j];
		if (multiplier_j != 0.) {
			for (std::map<int,double>::iterator it = current_iterate.constraints_jacobian[j].begin(); it != current_iterate.constraints_jacobian[j].end(); it++) {
				int variable_index = it->first;
				double derivative = it->second;
				rhs[variable_index] += multiplier_j*derivative;
			}
		}
	}
	
	int current_multiplier = 0;
	for (int i = 0; i < problem.number_variables; i++) {
		if (variable_status[i] == BOUNDED_LOWER || variable_status[i] == BOUNDED_BOTH_SIDES) {
			rhs[i] += this->multipliers[i];
			rhs[problem.number_variables + current_multiplier] = -current_iterate.x[i] + variable_lb[i] + this->mu/this->multipliers[current_multiplier];
			current_multiplier++;
		}
		if (variable_status[i] == BOUNDED_UPPER || variable_status[i] == BOUNDED_BOTH_SIDES) {
			rhs[i] += this->multipliers[i];
			rhs[problem.number_variables + current_multiplier] = -current_iterate.x[i] + variable_ub[i] + this->mu/this->multipliers[current_multiplier];
			current_multiplier++;
		}
	}
	
	for (int j = 0; j < problem.number_constraints; j++) {
		if (problem.constraint_status[j] == EQUAL_BOUNDS) {
			rhs[problem.number_variables + current_multiplier] = -current_iterate.constraints[j];
			current_multiplier++;
		}
		if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			rhs[problem.number_variables + current_multiplier] = -current_iterate.constraints[j] + problem.constraint_lb[j] + this->mu/this->multipliers[current_multiplier];
			current_multiplier++;
		}
		if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			rhs[problem.number_variables + current_multiplier] = -current_iterate.constraints[j] + problem.constraint_ub[j] + this->mu/this->multipliers[current_multiplier];
			current_multiplier++;
		}
	}
	return rhs;
}

std::vector<double> InteriorPoint::compute_slack_displacements(Problem& problem, double mu, std::vector<double>& solution) {
	std::vector<double> delta_s(slacks.size());
	int current_multiplier = 0;
	int current_slack = 0;
	for (int i = 0; i < problem.number_variables; i++) {
		if (this->variable_status[i] == BOUNDED_LOWER || this->variable_status[i] == BOUNDED_BOTH_SIDES) {
			delta_s[current_slack] = mu/this->multipliers[current_multiplier] - this->slacks[current_slack] - this->slacks[current_slack]/this->multipliers[current_multiplier]*solution[problem.number_variables + current_multiplier];
			current_multiplier++;
			current_slack++;
		}
		if (this->variable_status[i] == BOUNDED_UPPER || this->variable_status[i] == BOUNDED_BOTH_SIDES) {
			delta_s[current_slack] = -mu/this->multipliers[current_multiplier] - this->slacks[current_slack] - this->slacks[current_slack]/this->multipliers[current_multiplier]*solution[problem.number_variables + current_multiplier];
			current_multiplier++;
			current_slack++;
		}
	}
	for (int j = 0; j < problem.number_constraints; j++) {
		if (problem.constraint_status[j] == EQUAL_BOUNDS) {
			// no slack
			current_multiplier++;
		}
		if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			delta_s[current_slack] = mu/this->multipliers[current_multiplier] - this->slacks[current_slack] - this->slacks[current_slack]/this->multipliers[current_multiplier]*solution[problem.number_variables + current_multiplier];
			current_multiplier++;
			current_slack++;
		}
		if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			delta_s[current_slack] = -mu/this->multipliers[current_multiplier] - this->slacks[current_slack] - this->slacks[current_slack]/this->multipliers[current_multiplier]*solution[problem.number_variables + current_multiplier];
			current_multiplier++;
			current_slack++;
		}
	}
	return delta_s;
}
