#include "InteriorPoint.hpp"

InteriorPoint::InteriorPoint(): LocalApproximation("IPM") {
}

void InteriorPoint::allocate_solver(int number_variables, int number_constraints) {
}
		
LocalSolution InteriorPoint::compute_optimality_step(Problem& problem, Iterate& current_iterate, double objective_multiplier, double radius) {
	double tau = 0.995;
	
	/* initialize the multipliers and count the slacks */
	int number_variable_slacks = 0;
	int number_constraint_slacks = 0;
	std::vector<double> multipliers;
	
	double initial_multiplier = 1e-6;
	
	for (int i = 0; i < problem.number_variables; i++) {
		// EQUALITY, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES, UNBOUNDED
		if (problem.variable_status[i] == EQUAL_BOUNDS) {
			std::cout << "x" << i << ": EQUAL_BOUNDS\n";
		}
		else if (problem.variable_status[i] == BOUNDED_LOWER || problem.variable_status[i] == BOUNDED_BOTH_SIDES) {
			std::cout << "x" << i << ": BOUNDED_LOWER\n";
			number_variable_slacks++;
			multipliers.push_back(initial_multiplier); // positive multiplier
		}
		else if (problem.variable_status[i] == BOUNDED_UPPER || problem.variable_status[i] == BOUNDED_BOTH_SIDES) {
			std::cout << "x" << i << ": BOUNDED_UPPER\n";
			number_variable_slacks++;
			multipliers.push_back(-initial_multiplier); // negative multiplier
		}
	}
	for (int j = 0; j < problem.number_constraints; j++) {
		// EQUALITY, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES, UNBOUNDED
		if (problem.constraint_status[j] == EQUAL_BOUNDS) {
			std::cout << "c" << j << ": EQUALITY\n";
		}
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
	
	/* initialize the slacks */
	std::vector<double> s(number_slacks);
	for (int j = 0; j < number_slacks; j++) {
		s[j] = 1.;
	}
	
	std::vector<double> x = {0.48, 2.2};
	//std::vector<double> x = {0.5, 2};
	//std::vector<double> x = {-2, 1};
	s[0] = 1e-6;
	s[1] = 1e-6;
	s[2] = 4.5;
	multipliers[0] = -1000;
	multipliers[1] = 600.;
	multipliers[2] = 1e-6;
	double mu = 1e-9;

	std::cout << "\n";

	std::cout << "x is: (" << x[0] << ", " << x[1] << ")\n";
	std::cout << "s is: (" << s[0] << ", " << s[1] << ", " << s[2] << ")\n";
	std::cout << "λ is: (" << multipliers[0] << ", " << multipliers[1] << ", " << multipliers[2] << ")\n";
	std::cout << "mu is " << mu << "\n";

	
	/********************************/
	/* reduced primal-dual approach */
	/********************************/
	
	/* evaluate the Lagrangian Hessian */
	// recreate the multipliers
	std::vector<double> original_multipliers(problem.number_constraints);
	int current_slack = number_variable_slacks;
	for (int j = 0; j < problem.number_constraints; j++) {
		if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			original_multipliers[j] += multipliers[current_slack];
			current_slack++;
		}
		if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			original_multipliers[j] += multipliers[current_slack];
			current_slack++;
		}
	}
	// create the COO matrix
	CSCMatrix hessian = problem.lagrangian_hessian(x, objective_multiplier, original_multipliers);
	COOMatrix matrix = hessian.to_COO();
	matrix.size += number_slacks;
	
	/* fill the rest of the matrix */
	current_slack = problem.number_variables;
	// variables
	for (int i = 0; i < problem.number_variables; i++) {
		// EQUALITY, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES, UNBOUNDED
		if (problem.variable_status[i] == BOUNDED_LOWER || problem.variable_status[i] == BOUNDED_BOTH_SIDES) {
			matrix.add_term(1., i, current_slack);
			current_slack++;
		}
		if (problem.variable_status[i] == BOUNDED_UPPER || problem.variable_status[i] == BOUNDED_BOTH_SIDES) {
			matrix.add_term(1., i, current_slack);
			current_slack++;
		}
	}
		
	// fill with constraint jacobian
	// TODO duplicate if necessary
	std::vector<std::map<int,double> > constraints_sparse_jacobian = problem.constraints_sparse_jacobian(x);
	for (int j = 0; j < problem.number_constraints; j++) {
		for (std::map<int,double>::iterator it = constraints_sparse_jacobian[j].begin(); it != constraints_sparse_jacobian[j].end(); it++) {
			int variable_index = it->first;
			double derivative = it->second;
			matrix.add_term(derivative, variable_index, current_slack);
		}
		current_slack++;
	}
	
	/* diagonal terms -L^-1 S */
	for (int j = 0; j < number_slacks; j++) {
		matrix.add_term(-s[j]/multipliers[j], problem.number_variables + j, problem.number_variables + j);
	}

	/* right hand side */
	std::vector<double> rhs(problem.number_variables + number_slacks);
	// objective gradient
	std::map<int,double> objective_gradient = problem.objective_sparse_gradient(x);
	for (std::map<int,double>::iterator it = objective_gradient.begin(); it != objective_gradient.end(); it++) {
		int variable_index = it->first;
		double derivative = it->second;
		/* scale the objective gradient */
		rhs[variable_index] -= objective_multiplier*derivative;
	}
	// constraints
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
	
	current_slack = 0;
	for (int i = 0; i < problem.number_variables; i++) {
		// EQUALITY, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES, UNBOUNDED
		if (problem.variable_status[i] == BOUNDED_LOWER || problem.variable_status[i] == BOUNDED_BOTH_SIDES) {
			rhs[i] += multipliers[i];
			rhs[problem.number_variables + current_slack] = -x[i] + problem.variable_lb[i] + mu/multipliers[current_slack];
			current_slack++;
		}
		if (problem.variable_status[i] == BOUNDED_UPPER || problem.variable_status[i] == BOUNDED_BOTH_SIDES) {
			rhs[i] += multipliers[i];
			rhs[problem.number_variables + current_slack] = -x[i] + problem.variable_ub[i] - mu/multipliers[current_slack] - 2*s[current_slack];
			current_slack++;
		}
	}
	for (int j = 0; j < problem.number_constraints; j++) {
		if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			rhs[problem.number_variables + current_slack] = -problem.evaluate_constraint(j, x) + problem.constraint_lb[j] + mu/multipliers[current_slack];
			current_slack++;
		}
		if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			rhs[problem.number_variables + current_slack] = -problem.evaluate_constraint(j, x) + problem.constraint_ub[j] - mu/multipliers[current_slack] - 2*s[current_slack];
			current_slack++;
		}
	}
	
	LocalSolution solution = this->solver.solve(matrix, rhs);
	
	/* negate multipliers displacements (Nocedal p590) */
	for (int j = problem.number_variables; j < problem.number_variables + number_slacks; j++) {
		solution.x[j] = -solution.x[j];
	}
	
	std::cout << "solution:\n";
	for (int k = 0; k < problem.number_variables; k++) {
		std::cout << "Δx" << k << " = " << solution.x[k] << "\n";
	}
	for (int k = 0; k < number_slacks; k++) {
		std::cout << "Δλ" << k << " = " << solution.x[problem.number_variables+k] << "\n";
	}
	/* recreate slack displacements */
	std::vector<double> new_delta_s(number_slacks);
	current_slack = 0;
	for (int i = 0; i < problem.number_variables; i++) {
		// EQUALITY, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES, UNBOUNDED
		if (problem.variable_status[i] == BOUNDED_LOWER || problem.variable_status[i] == BOUNDED_BOTH_SIDES) {
			new_delta_s[current_slack] = mu/multipliers[current_slack] - s[current_slack] - s[current_slack]/multipliers[current_slack]*solution.x[problem.number_variables + current_slack];
			current_slack++;
		}
		if (problem.variable_status[i] == BOUNDED_UPPER || problem.variable_status[i] == BOUNDED_BOTH_SIDES) {
			new_delta_s[current_slack] = mu/multipliers[current_slack] + s[current_slack] + s[current_slack]/multipliers[current_slack]*solution.x[problem.number_variables + current_slack];
			current_slack++;
		}
	}
	for (int j = 0; j < problem.number_constraints; j++) {
		// EQUALITY, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES, UNBOUNDED
		if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			new_delta_s[current_slack] = mu/multipliers[current_slack] - s[current_slack] - s[current_slack]/multipliers[current_slack]*solution.x[problem.number_variables + current_slack];
			current_slack++;
		}
		if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			new_delta_s[current_slack] = mu/multipliers[current_slack] + s[current_slack] + s[current_slack]/multipliers[current_slack]*solution.x[problem.number_variables + current_slack];
			current_slack++;
		}
	}

	for (int j = 0; j < number_slacks; j++) {
		std::cout << "Δs" << j << " = " << new_delta_s[j] << "\n";
	}
	std::cout << "\n";
	
	/* fraction to boundary rule */
	double alpha_x = 1.;
	for (int j = 0; j < number_slacks; j++) {
		/* sj */
		double trial_alpha_sj = -tau*s[j]/new_delta_s[j];
		if (0 < trial_alpha_sj && trial_alpha_sj <= 1.) {
			alpha_x = std::min(alpha_x, trial_alpha_sj);
		}
	}
	std::cout << "alpha_x = " << alpha_x << "\n";
	
	double alpha_mult = 1.;
	for (int j = 0; j < number_slacks; j++) {
		/* lambda_j */
		double trial_alpha_lj = -tau*multipliers[j]/solution.x[problem.number_variables+j];
		if (0 < trial_alpha_lj && trial_alpha_lj <= 1.) {
			alpha_mult = std::min(alpha_mult, trial_alpha_lj);
		}
	}
	std::cout << "alpha_λ = " << alpha_mult << "\n\n";

	std::cout << "New x: (";
	for (int i = 0; i < problem.number_variables; i++) {
		x[i] += alpha_x*solution.x[i];
		std::cout << x[i] << ", ";
	}
	std::cout << ")\n";
	std::cout << "New s: (";
	for (int j = 0; j < number_slacks; j++) {
		s[j] += alpha_x*new_delta_s[j];
		std::cout << s[j] << ", ";
	}
	std::cout << ")\n";
	std::cout << "New λ: (";
	for (int j = 0; j < number_slacks; j++) {
		multipliers[j] += alpha_mult*solution.x[problem.number_variables+j];
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
