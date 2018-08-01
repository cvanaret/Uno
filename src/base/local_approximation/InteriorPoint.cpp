#include "InteriorPoint.hpp"

InteriorPoint::InteriorPoint(): LocalApproximation("IPM") {
}

void InteriorPoint::allocate_solver(int number_variables, int number_constraints) {
}
		
LocalSolution InteriorPoint::compute_optimality_step(Problem& problem, Iterate& current_iterate, double objective_multiplier, double radius) {
	/* reduced primal-dual approach */
	double tau = 0.995;
	
	/* initialize the multipliers and count the slacks */
	int number_variable_slacks = 0;
	int number_constraint_slacks = 0;
	std::vector<double> multipliers;
	
	double initial_multiplier = 1e-6;
	
	/* TODO: handle the trust region */
	
	for (int i = 0; i < problem.number_variables; i++) {
		if (problem.variable_status[i] == EQUAL_BOUNDS) {
			std::cout << "x" << i << ": EQUAL_BOUNDS\n";
			throw std::out_of_range("IPM: variable has equal bounds, behavior not specified");
		}
		if (problem.variable_status[i] == BOUNDED_LOWER || problem.variable_status[i] == BOUNDED_BOTH_SIDES) {
			std::cout << "x" << i << ": BOUNDED_LOWER\n";
			number_variable_slacks++;
			multipliers.push_back(initial_multiplier); // positive multiplier
		}
		if (problem.variable_status[i] == BOUNDED_UPPER || problem.variable_status[i] == BOUNDED_BOTH_SIDES) {
			std::cout << "x" << i << ": BOUNDED_UPPER\n";
			number_variable_slacks++;
			multipliers.push_back(-initial_multiplier); // negative multiplier
		}
	}
	for (int j = 0; j < problem.number_constraints; j++) {
		if (problem.constraint_status[j] == EQUAL_BOUNDS) {
			std::cout << "c" << j << ": EQUAL_BOUNDS\n";
			// no slack
			multipliers.push_back(initial_multiplier); // unconstrained multiplier
		}
		/* slacks for inequality constraints */
		if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			std::cout << "c" << j << ": BOUNDED_LOWER\n";
			number_constraint_slacks++;
			multipliers.push_back(initial_multiplier); // positive multiplier
		}
		if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			std::cout << "c" << j << ": BOUNDED_UPPER\n";
			number_constraint_slacks++;
			multipliers.push_back(-initial_multiplier); // negative multiplier
		}
	}
	int number_slacks = number_variable_slacks + number_constraint_slacks;
	std::cout << number_slacks << " slack variables\n";
	std::cout << multipliers.size() << " constraints/multipliers\n";
	
	/* initialize the slacks */
	std::vector<double> s(number_slacks);
	for (int j = 0; j < number_slacks; j++) {
		s[j] = 1.;
	}
	
	std::vector<double> x(current_iterate.x);
	//std::vector<double> x = {0.48, 2.2};
	////std::vector<double> x = {0.5, 2};
	////std::vector<double> x = {-2, 1};
	//s[0] = 1e-6;
	//s[1] = 1e-6;
	//s[2] = 4.5;
	//multipliers[0] = -1000;
	//multipliers[1] = 600.;
	//multipliers[2] = 1e-6;
	double mu = 10.;

	std::cout << "\n";

	std::cout << "x is: "; print_vector(std::cout, x);
	std::cout << "s is: "; print_vector(std::cout, s);
	std::cout << "λ is: "; print_vector(std::cout, multipliers);
	std::cout << "mu is " << mu << "\n\n";

	
	/***************************/
	/* sparse symmetric matrix */
	/***************************/
	
	/* evaluate the Lagrangian Hessian */
	// recreate the multipliers
	std::vector<double> original_multipliers(problem.number_constraints);
	int current_multiplier = number_variable_slacks;
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
	COOMatrix matrix = hessian.to_COO();
	matrix.size += multipliers.size();
	
	/* Jacobian of bound constraints */
	int current_column = problem.number_variables;
	for (int i = 0; i < problem.number_variables; i++) {
		if (problem.variable_status[i] != UNBOUNDED) {
			matrix.add_term(1., i, current_column);
			current_column++;
		}
		if (problem.variable_status[i] == BOUNDED_BOTH_SIDES) {
			matrix.add_term(1., i, current_column);
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
				matrix.add_term(derivative, variable_index, current_column);
			}
			current_column++;
		}
	}
	
	/* variable diagonal terms -L^-1 S */
	int current_slack = 0;
	current_multiplier = 0;
	for (int i = 0; i < problem.number_variables; i++) {
		if (problem.variable_status[i] != UNBOUNDED) { 
			current_column = problem.number_variables + current_multiplier;
			matrix.add_term(-s[current_slack]/multipliers[current_multiplier], current_column, current_column);
			current_slack++;
			current_multiplier++;
		}
		if (problem.variable_status[i] == BOUNDED_BOTH_SIDES) {
			current_column = problem.number_variables + current_multiplier;
			matrix.add_term(-s[current_slack]/multipliers[current_multiplier], current_column, current_column);
			current_slack++;
			current_multiplier++;
		}
	}
	/* constraint diagonal terms -L^-1 S */
	for (int j = 0; j < problem.number_constraints; j++) {
		if (problem.constraint_status[j] == EQUAL_BOUNDS) {
			// no slack
			current_multiplier++;
		}
		else {
			current_column = problem.number_variables + current_multiplier;
			matrix.add_term(-s[current_slack]/multipliers[current_multiplier], current_column, current_column);
			current_slack++;
			current_multiplier++;
			if (problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
				current_column = problem.number_variables + current_multiplier;
				matrix.add_term(-s[current_slack]/multipliers[current_multiplier], current_column, current_column);
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
		/* scale the objective gradient */
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
	current_slack = 0;
	for (int i = 0; i < problem.number_variables; i++) {
		if (problem.variable_status[i] == BOUNDED_LOWER || problem.variable_status[i] == BOUNDED_BOTH_SIDES) {
			rhs[i] += multipliers[i];
			rhs[problem.number_variables + current_multiplier] = -x[i] + problem.variable_lb[i] + mu/multipliers[current_multiplier];
			current_multiplier++;
			current_slack++;
		}
		if (problem.variable_status[i] == BOUNDED_UPPER || problem.variable_status[i] == BOUNDED_BOTH_SIDES) {
			rhs[i] += multipliers[i];
			rhs[problem.number_variables + current_multiplier] = -x[i] + problem.variable_ub[i] - mu/multipliers[current_multiplier] - 2*s[current_slack];
			current_multiplier++;
			current_slack++;
		}
	}
	
	for (int j = 0; j < problem.number_constraints; j++) {
		if (problem.constraint_status[j] == EQUAL_BOUNDS) {
			rhs[problem.number_variables + current_multiplier] = -problem.evaluate_constraint(j, x);
			current_multiplier++;
			// no slack
		}
		if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			rhs[problem.number_variables + current_multiplier] = -problem.evaluate_constraint(j, x) + problem.constraint_lb[j] + mu/multipliers[current_multiplier];
			current_multiplier++;
			current_slack++;
		}
		if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			rhs[problem.number_variables + current_multiplier] = -problem.evaluate_constraint(j, x) + problem.constraint_ub[j] - mu/multipliers[current_multiplier] - 2*s[current_slack];
			current_multiplier++;
			current_slack++;
		}
	}
	
	/************/
	/* solution */
	/************/
	
	/* compute the solution (Δx, -Δλ) */
	LocalSolution solution = this->solver.solve(matrix, rhs);
	
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
	std::vector<double> delta_s(number_slacks);
	current_multiplier = 0;
	current_slack = 0;
	for (int i = 0; i < problem.number_variables; i++) {
		if (problem.variable_status[i] == BOUNDED_LOWER || problem.variable_status[i] == BOUNDED_BOTH_SIDES) {
			delta_s[current_slack] = mu/multipliers[current_multiplier] - s[current_slack] - s[current_slack]/multipliers[current_multiplier]*solution.x[problem.number_variables + current_multiplier];
			current_multiplier++;
			current_slack++;
		}
		if (problem.variable_status[i] == BOUNDED_UPPER || problem.variable_status[i] == BOUNDED_BOTH_SIDES) {
			delta_s[current_slack] = mu/multipliers[current_multiplier] + s[current_slack] + s[current_slack]/multipliers[current_multiplier]*solution.x[problem.number_variables + current_multiplier];
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
			delta_s[current_slack] = mu/multipliers[current_multiplier] - s[current_slack] - s[current_slack]/multipliers[current_multiplier]*solution.x[problem.number_variables + current_multiplier];
			current_multiplier++;
			current_slack++;
		}
		if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			delta_s[current_slack] = mu/multipliers[current_multiplier] + s[current_slack] + s[current_slack]/multipliers[current_multiplier]*solution.x[problem.number_variables + current_multiplier];
			current_multiplier++;
			current_slack++;
		}
	}

	for (int j = 0; j < number_slacks; j++) {
		std::cout << "Δs" << j << " = " << delta_s[j] << "\n";
	}
	std::cout << "\n";
	
	/* fraction to boundary rule for slacks */
	double alpha_s = 1.;
	for (int j = 0; j < number_slacks; j++) {
		double trial_alpha_sj = -tau*s[j]/delta_s[j];
		if (0 < trial_alpha_sj && trial_alpha_sj <= 1.) {
			alpha_s = std::min(alpha_s, trial_alpha_sj);
		}
	}
	std::cout << "alpha_s = " << alpha_s << "\n";
	
	/* fraction to boundary rule for multipliers */
	double alpha_mult = 1.;
	current_multiplier = 0;
	for (int j = 0; j < problem.number_constraints; j++) {
		if (problem.constraint_status[j] == EQUAL_BOUNDS) {
			// multiplier is unconstrained
			current_multiplier++;
		}
		else {
			int number_copies = (problem.variable_status[j] == BOUNDED_BOTH_SIDES) ? 2 : 1;
			for (int copy = 1; copy <= number_copies; copy++) {
				double trial_alpha_lj = -tau*multipliers[current_multiplier]/solution.x[problem.number_variables + current_multiplier];
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
	for (int j = 0; j < number_slacks; j++) {
		s[j] += alpha_s*delta_s[j];
		std::cout << s[j] << ", ";
	}
	std::cout << ")\n";
	std::cout << "New λ: (";
	for (unsigned int j = 0; j < multipliers.size(); j++) {
		multipliers[j] += alpha_mult*solution.x[problem.number_variables + j];
		std::cout << multipliers[j] << ", ";
	}
	std::cout << ")\n";
	
	return solution;
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
