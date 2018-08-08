#include "InteriorPoint.hpp"

InteriorPoint::InteriorPoint(): Subproblem("IPM"), tau(0.995), default_multiplier(0.5), default_slack(1.) {
}

void InteriorPoint::initialize(Problem& problem, Iterate& current_iterate, int number_variables, int number_constraints, bool use_trust_region) {
}
		
LocalSolution InteriorPoint::compute_optimality_step(Problem& problem, Iterate& current_iterate, double objective_multiplier, double radius) {
	/* reduced primal-dual approach */
	std::vector<double> x(current_iterate.x);
	
	/* initialize the multipliers and count the slacks */
	std::vector<double> multipliers;
	std::vector<double> slacks;

	/* build trust region */
	std::vector<double> variable_lb(problem.number_variables);
	std::vector<double> variable_ub(problem.number_variables);
	for (int i = 0; i < problem.number_variables; i++) {
		variable_lb[i] = std::max(x[i] - radius, problem.variable_lb[i]);
		variable_ub[i] = std::min(x[i] + radius, problem.variable_ub[i]);
	}
	for (int i = 0; i < problem.number_variables; i++) {
		std::cout << "x" << i << " in [" << variable_lb[i] << ", " << variable_ub[i] << "]\n";
	}
	std::vector<ConstraintType> variable_status = problem.determine_constraints_types(variable_lb, variable_ub);
	
	/* determine the slacks and multipliers */
	for (int i = 0; i < problem.number_variables; i++) {
		if (variable_status[i] == EQUAL_BOUNDS) {
			std::cout << "x" << i << ": EQUAL_BOUNDS\n";
			throw std::out_of_range("IPM: variable has equal bounds, behavior not specified");
		}
		if (variable_status[i] == BOUNDED_LOWER || variable_status[i] == BOUNDED_BOTH_SIDES) {
			std::cout << "x" << i << ": BOUNDED_LOWER\n";
			multipliers.push_back(this->default_multiplier); // positive multiplier
			slacks.push_back(x[i] - variable_lb[i]); // exact slack
		}
		if (variable_status[i] == BOUNDED_UPPER || variable_status[i] == BOUNDED_BOTH_SIDES) {
			std::cout << "x" << i << ": BOUNDED_UPPER\n";
			multipliers.push_back(-this->default_multiplier); // negative multiplier
			slacks.push_back(variable_ub[i] - x[i]); // exact slack
		}
	}
	int number_variable_multipliers = multipliers.size();
	
	for (int j = 0; j < problem.number_constraints; j++) {
		if (problem.constraint_status[j] == EQUAL_BOUNDS) {
			std::cout << "c" << j << ": EQUAL_BOUNDS\n";
			multipliers.push_back(this->default_multiplier); // unconstrained multiplier
			// no slack
		}
		/* slacks for inequality constraints */
		if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			std::cout << "c" << j << ": BOUNDED_LOWER\n";
			multipliers.push_back(this->default_multiplier); // positive multiplier
			slacks.push_back(this->default_slack);
		}
		if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			std::cout << "c" << j << ": BOUNDED_UPPER\n";
			multipliers.push_back(-this->default_multiplier); // negative multiplier
			slacks.push_back(this->default_slack);
		}
	}
	
	std::cout << slacks.size() << " slack variables\n";
	std::cout << multipliers.size() << " constraints/multipliers\n";
	
	//x = {0.48, 2.2};
	//x = {0.5, 2};
	//std::vector<double> x = {-2, 1};
	//s[0] = 1e-6;
	//s[1] = 1e-6;
	//s[2] = 4.5;
	//multipliers[0] = -1000;
	//multipliers[1] = 600.;
	//multipliers[2] = 1e-6;
	double mu = 1e-6;
	
	std::cout << "\n";

	std::cout << "x is: "; print_vector(std::cout, x);
	std::cout << "s is: "; print_vector(std::cout, slacks);
	std::cout << "λ is: "; print_vector(std::cout, multipliers);
	std::cout << "mu is " << mu << "\n\n";

	//for (int iter = 1; iter < 3; iter++) {
	
	/***************************/
	/* sparse symmetric matrix */
	/***************************/
	
	/* evaluate the Lagrangian Hessian */
	// recreate the multipliers
	std::vector<double> original_multipliers(problem.number_constraints);
	int current_multiplier = number_variable_multipliers;
	for (int j = 0; j < problem.number_constraints; j++) {
		original_multipliers[j] = multipliers[current_multiplier];
		current_multiplier++;
		// if two bounds, add the two multipliers
		if (problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			original_multipliers[j] += multipliers[current_multiplier];
			current_multiplier++;
		}
	}
	// create the COO Hessian matrix
	CSCMatrix hessian = problem.lagrangian_hessian(x, objective_multiplier, original_multipliers);
	COOMatrix kkt_matrix = hessian.to_COO();
	kkt_matrix.size += multipliers.size();
	int current_column = problem.number_variables;
	
	/* Jacobian of bound constraints */
	for (int i = 0; i < problem.number_variables; i++) {
		if (variable_status[i] != UNBOUNDED) {
			kkt_matrix.add_term(1., i, current_column);
			current_column++;
		}
		if (variable_status[i] == BOUNDED_BOTH_SIDES) {
			kkt_matrix.add_term(1., i, current_column);
			current_column++;
		}
	}
	/* Jacobian of general constraints */
	std::vector<std::map<int,double> > constraints_sparse_jacobian = problem.constraints_sparse_jacobian(x);
	for (int j = 0; j < problem.number_constraints; j++) {
		int number_copies = (problem.constraint_status[j] == BOUNDED_BOTH_SIDES) ? 2 : 1;
		for (int copy = 1; copy <= number_copies; copy++) {
			for (std::map<int,double>::iterator it = constraints_sparse_jacobian[j].begin(); it != constraints_sparse_jacobian[j].end(); it++) {
				int variable_index = it->first;
				double derivative = it->second;
				kkt_matrix.add_term(derivative, variable_index, current_column);
			}
			current_column++;
		}
	}
	
	/* variable diagonal terms +/-L^-1 S */
	int current_slack = 0;
	current_multiplier = 0;
	for (int i = 0; i < problem.number_variables; i++) {
		if (variable_status[i] == BOUNDED_LOWER || variable_status[i] == BOUNDED_BOTH_SIDES) {
			current_column = problem.number_variables + current_multiplier;
			kkt_matrix.add_term(-slacks[current_slack]/multipliers[current_multiplier], current_column, current_column);
			current_slack++;
			current_multiplier++;
		}
		if (variable_status[i] == BOUNDED_UPPER || variable_status[i] == BOUNDED_BOTH_SIDES) {
			current_column = problem.number_variables + current_multiplier;
			kkt_matrix.add_term(slacks[current_slack]/multipliers[current_multiplier], current_column, current_column);
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
				kkt_matrix.add_term(-slacks[current_slack]/multipliers[current_multiplier], current_column, current_column);
				current_slack++;
				current_multiplier++;
			}
			if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
				current_column = problem.number_variables + current_multiplier;
				kkt_matrix.add_term(slacks[current_slack]/multipliers[current_multiplier], current_column, current_column);
				current_slack++;
				current_multiplier++;
			}
		}
	}
	
	/*******************/
	/* right-hand side */
	/*******************/
	
	std::vector<double> rhs(problem.number_variables + multipliers.size());
	
	/* objective gradient */
	std::map<int,double> objective_gradient = problem.objective_sparse_gradient(x);
	for (std::map<int,double>::iterator it = objective_gradient.begin(); it != objective_gradient.end(); it++) {
		int variable_index = it->first;
		double derivative = it->second;
		rhs[variable_index] -= objective_multiplier*derivative;
	}
	/* constraint gradients */
	for (int j = 0; j < problem.number_constraints; j++) {
		double multiplier_j = original_multipliers[j];
		if (multiplier_j != 0.) {
			for (std::map<int,double>::iterator it = constraints_sparse_jacobian[j].begin(); it != constraints_sparse_jacobian[j].end(); it++) {
				int variable_index = it->first;
				double derivative = it->second;
				rhs[variable_index] += multiplier_j*derivative;
			}
		}
	}
	
	current_multiplier = 0;
	for (int i = 0; i < problem.number_variables; i++) {
		if (variable_status[i] == BOUNDED_LOWER || variable_status[i] == BOUNDED_BOTH_SIDES) {
			rhs[i] += multipliers[i];
			rhs[problem.number_variables + current_multiplier] = -x[i] + variable_lb[i] + mu/multipliers[current_multiplier];
			current_multiplier++;
		}
		if (variable_status[i] == BOUNDED_UPPER || variable_status[i] == BOUNDED_BOTH_SIDES) {
			rhs[i] += multipliers[i];
			rhs[problem.number_variables + current_multiplier] = -x[i] + variable_ub[i] + mu/multipliers[current_multiplier];
			std::cout << "rhs[" << problem.number_variables + current_multiplier << "] = -x[" << i << "] + " << variable_ub[i] << " + " << mu << "/" << multipliers[current_multiplier] << "\n";
			current_multiplier++;
		}
	}
	
	for (int j = 0; j < problem.number_constraints; j++) {
		if (problem.constraint_status[j] == EQUAL_BOUNDS) {
			rhs[problem.number_variables + current_multiplier] = -current_iterate.constraints[j];
			current_multiplier++;
		}
		if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			rhs[problem.number_variables + current_multiplier] = -current_iterate.constraints[j] + problem.constraint_lb[j] + mu/multipliers[current_multiplier];
			current_multiplier++;
		}
		if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			rhs[problem.number_variables + current_multiplier] = -current_iterate.constraints[j] + problem.constraint_ub[j] + mu/multipliers[current_multiplier];
			std::cout << "rhs[" << problem.number_variables + current_multiplier << "] = " << -current_iterate.constraints[j] << " + " << problem.constraint_ub[j] << " + " << mu << "/" << multipliers[current_multiplier] << "\n";
			current_multiplier++;
		}
	}
	
	for (int k = 0; k < kkt_matrix.matrix.size(); k++) {
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
	std::vector<double> delta_s = this->compute_slack_displacements(problem, variable_status, mu, multipliers, slacks, solution.x);

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
	for (int i = 0; i < problem.number_variables; i++) {
		x[i] += alpha_s*solution.x[i];
		std::cout << x[i] << ", ";
	}
	std::cout << ")\n";
	std::cout << "New s: (";
	for (unsigned int j = 0; j < slacks.size(); j++) {
		slacks[j] += alpha_s*delta_s[j];
		std::cout << slacks[j] << ", ";
	}
	std::cout << ")\n";
	std::cout << "New λ: (";
	for (unsigned int j = 0; j < multipliers.size(); j++) {
		multipliers[j] += alpha_mult*solution.x[problem.number_variables + j];
		std::cout << multipliers[j] << ", ";
	}
	std::cout << ")\n";
	
	/*******************************/
	/* update of barrier parameter */
	/*******************************/
	
	double sum = 0.;
	current_multiplier = 0;
	current_slack = 0;
	
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
	mu = sum;
	//}
	
	LocalSolution sol(x, 0, 0);
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

std::vector<double> InteriorPoint::compute_slack_displacements(Problem& problem, std::vector<ConstraintType>& variable_status, double mu, std::vector<double>& multipliers, std::vector<double>& slacks, std::vector<double>& solution) {
	std::vector<double> delta_s(slacks.size());
	int current_multiplier = 0;
	int current_slack = 0;
	for (int i = 0; i < problem.number_variables; i++) {
		if (variable_status[i] == BOUNDED_LOWER || variable_status[i] == BOUNDED_BOTH_SIDES) {
			delta_s[current_slack] = mu/multipliers[current_multiplier] - slacks[current_slack] - slacks[current_slack]/multipliers[current_multiplier]*solution[problem.number_variables + current_multiplier];
			current_multiplier++;
			current_slack++;
		}
		if (variable_status[i] == BOUNDED_UPPER || variable_status[i] == BOUNDED_BOTH_SIDES) {
			delta_s[current_slack] = -mu/multipliers[current_multiplier] - slacks[current_slack] - slacks[current_slack]/multipliers[current_multiplier]*solution[problem.number_variables + current_multiplier];
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
			delta_s[current_slack] = mu/multipliers[current_multiplier] - slacks[current_slack] - slacks[current_slack]/multipliers[current_multiplier]*solution[problem.number_variables + current_multiplier];
			current_multiplier++;
			current_slack++;
		}
		if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			delta_s[current_slack] = -mu/multipliers[current_multiplier] - slacks[current_slack] - slacks[current_slack]/multipliers[current_multiplier]*solution[problem.number_variables + current_multiplier];
			current_multiplier++;
			current_slack++;
		}
	}
	return delta_s;
}
