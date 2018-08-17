#include <cmath>
#include "InteriorPoint.hpp"

InteriorPoint::InteriorPoint(): Subproblem("IPM"), mu(0.1), number_slacks(0), tau_min(0.99),
		default_multiplier(1.), k_sigma(1e10), inertia_term(0.) {
}

void InteriorPoint::initialize(Problem& problem, Iterate& current_iterate, int number_variables, int number_constraints, bool use_trust_region) {
	/* if trust region is used, bound constraints become range constraints */
	this->variable_status = problem.variable_status;
	if (use_trust_region) {
		for (int i = 0; i < problem.number_variables; i++) {
			if (this->variable_status[i] != EQUAL_BOUNDS) {
				this->variable_status[i] = BOUNDED_BOTH_SIDES;
			}
		}
	}
	
	/* determine the variable multipliers */
	for (int i = 0; i < problem.number_variables; i++) {
		if (this->variable_status[i] == BOUNDED_LOWER || this->variable_status[i] == BOUNDED_BOTH_SIDES) {
			this->bound_multipliers.push_back(this->default_multiplier); // positive multiplier
		}
		if (this->variable_status[i] == BOUNDED_UPPER || this->variable_status[i] == BOUNDED_BOTH_SIDES) {
			this->bound_multipliers.push_back(-this->default_multiplier); // negative multiplier
		}
	}
	/* determine the inequality constraint slacks */
	for (int j = 0; j < problem.number_constraints; j++) {
		if (problem.constraint_status[j] != EQUAL_BOUNDS) {
			double slack_value = this->project_variable_in_bounds(current_iterate.constraints[j], problem.constraint_lb[j], problem.constraint_ub[j]);
			current_iterate.x.push_back(slack_value);
		}
	}
	
	this->number_slacks = current_iterate.x.size() - problem.number_variables;
	std::cout << this->number_slacks << " slacks\n";
	
	/* compute first-order information */
	current_iterate.compute_objective_gradient(problem);
	current_iterate.compute_constraints_jacobian(problem);
	/* compute least-square multipliers */
	this->constraint_multipliers = this->estimate_initial_multipliers(problem, current_iterate);
	
	std::cout << this->bound_multipliers.size() << " bound multipliers\n";
	std::cout << this->constraint_multipliers.size() << " constraint multipliers\n";
	return;
}

/* reduced primal-dual approach */
LocalSolution InteriorPoint::compute_optimality_step(Problem& problem, Iterate& current_iterate, double radius) {
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
	for (int j = 0; j < problem.number_constraints; j++) {
		std::cout << "c" << j << " in [" << problem.constraint_lb[j] << ", " << problem.constraint_ub[j] << "]\n";
	}
	std::cout << "\n";

	std::cout << "x/s is: "; print_vector(std::cout, current_iterate.x);
	std::cout << "λ is: "; print_vector(std::cout, this->constraint_multipliers);
	std::cout << "z is: "; print_vector(std::cout, this->bound_multipliers);
	std::cout << "mu is " << this->mu << "\n";
	std::cout << "Constraints: "; print_vector(std::cout, current_iterate.constraints);
	std::cout << "\n";
	
	/* compute first-order information */
	current_iterate.compute_objective_gradient(problem);
	current_iterate.compute_constraints_jacobian(problem);
	
	for (std::map<int,double>::iterator it = current_iterate.objective_gradient.begin(); it != current_iterate.objective_gradient.end(); it++) {
		int variable_index = it->first;
		double derivative = it->second;
		std::cout << "df/dx" << variable_index << " = " << derivative << "\n";
	}
	
	/*******************************/
	/* sparse symmetric KKT matrix */
	/*******************************/
	/* evaluate the Lagrangian Hessian */
	// recreate the original multipliers
	std::vector<double> original_multipliers(problem.number_constraints);
	int current_multiplier = 0;
	for (int j = 0; j < problem.number_constraints; j++) {
		original_multipliers[j] = this->constraint_multipliers[current_multiplier];
		current_multiplier++;
		// if both bounds, add the two multipliers
		if (problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			original_multipliers[j] += this->constraint_multipliers[current_multiplier];
			current_multiplier++;
		}
	}
	// create the KKT matrix
	COOMatrix kkt_matrix = this->generate_kkt_matrix(problem, current_iterate, original_multipliers, variable_lb, variable_ub);
	
	/*******************/
	/* right-hand side */
	/*******************/
	std::vector<double> rhs = this->generate_rhs(problem, current_iterate, original_multipliers, variable_lb, variable_ub);
	
	/************/
	/* solution */
	/************/
	/* compute the solution (Δx, -Δλ) */
	LocalSolution solution = this->solver.solve(kkt_matrix, rhs, this->data);
	
	/* retrieve +Δλ (Nocedal p590) */
	for (unsigned int j = problem.number_variables; j < problem.number_variables + constraint_multipliers.size(); j++) {
		solution.x[j] = -solution.x[j];
	}
	
	/* compute slack displacements Δs */
	std::vector<double> delta_s = this->compute_slack_displacements(problem, current_iterate, solution.x);

	/* compute bound multiplier displacements Δz */
	std::vector<double> delta_z = this->compute_bound_multiplier_displacements(problem, current_iterate, solution.x, variable_lb, variable_ub);

	std::cout << "MA57 solution:\n";
	for (int k = 0; k < problem.number_variables; k++) {
		std::cout << "Δx" << k << " = " << solution.x[k] << "\n";
	}
	for (unsigned int j = 0; j < this->constraint_multipliers.size(); j++) {
		std::cout << "Δλ" << j << " = " << solution.x[problem.number_variables + j] << "\n";
	}
	for (unsigned int j = 0; j < delta_s.size(); j++) {
		std::cout << "Δs" << j << " = " << delta_s[j] << "\n";
	}
	for (unsigned int j = 0; j < delta_z.size(); j++) {
		std::cout << "Δz" << j << " = " << delta_z[j] << "\n";
	}
	std::cout << "\n";
	
	/* fraction to boundary rule for variables */
	double tau = std::max(this->tau_min, 1-this->mu);
	double primal_length = 1.;
	for (int i = 0; i < problem.number_variables; i++) {
		if (variable_status[i] == BOUNDED_LOWER || variable_status[i] == BOUNDED_BOTH_SIDES) {
			double trial_alpha_xi = -tau*(current_iterate.x[i] - variable_lb[i])/solution.x[i];
			if (0 < trial_alpha_xi && trial_alpha_xi <= 1.) {
				primal_length = std::min(primal_length, trial_alpha_xi);
			}
		}
		if (variable_status[i] == BOUNDED_UPPER || variable_status[i] == BOUNDED_BOTH_SIDES) {
			double trial_alpha_xi = -tau*(current_iterate.x[i] - variable_ub[i])/solution.x[i];
			if (0 < trial_alpha_xi && trial_alpha_xi <= 1.) {
				primal_length = std::min(primal_length, trial_alpha_xi);
			}
		}
	}
	int current_slack = 0;
	for (int j = 0; j < problem.number_constraints; j++) {
		if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			double trial_primal_lengthj = -tau*(current_iterate.x[problem.number_variables + current_slack] - problem.constraint_lb[j])/delta_s[current_slack];
			if (0 < trial_primal_lengthj && trial_primal_lengthj <= 1.) {
				primal_length = std::min(primal_length, trial_primal_lengthj);
			}
		}
		if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			double trial_primal_lengthj = -tau*(current_iterate.x[problem.number_variables + current_slack] - problem.constraint_ub[j])/delta_s[current_slack];
			if (0 < trial_primal_lengthj && trial_primal_lengthj <= 1.) {
				primal_length = std::min(primal_length, trial_primal_lengthj);
			}
		}
		if (problem.constraint_status[j] != EQUAL_BOUNDS) {
			current_slack++;
		}
	}
	std::cout << "primal length = " << primal_length << "\n";
	
	/* fraction to boundary rule for multipliers */
	double dual_length = 1.;
	current_multiplier = 0;
	for (int i = 0; i < problem.number_variables; i++) {
		if (variable_status[i] == BOUNDED_LOWER || variable_status[i] == BOUNDED_BOTH_SIDES) {
			double trial_alpha_zj = -tau*this->bound_multipliers[current_multiplier]/delta_z[current_multiplier];
			if (0 < trial_alpha_zj && trial_alpha_zj <= 1.) {
				dual_length = std::min(dual_length, trial_alpha_zj);
			}
			current_multiplier++;
		}
		if (variable_status[i] == BOUNDED_UPPER || variable_status[i] == BOUNDED_BOTH_SIDES) {
			double trial_alpha_zj = -tau*this->bound_multipliers[current_multiplier]/delta_z[current_multiplier];
			if (0 < trial_alpha_zj && trial_alpha_zj <= 1.) {
				dual_length = std::min(dual_length, trial_alpha_zj);
			}
			current_multiplier++;
		}
	}
	std::cout << "dual length = " << dual_length << "\n\n";

	/* update */
	std::vector<double> trial_primal(current_iterate.x.size());
	for (int i = 0; i < problem.number_variables; i++) {
		trial_primal[i] = current_iterate.x[i] + primal_length*solution.x[i];
	}
	for (int j = 0; j < this->number_slacks; j++) {
		trial_primal[problem.number_variables + j] = current_iterate.x[problem.number_variables + j] + primal_length*delta_s[j];
	}
	std::cout << "New x/s: "; print_vector(std::cout, trial_primal);
	
	for (unsigned int j = 0; j < this->constraint_multipliers.size(); j++) {
		this->constraint_multipliers[j] += primal_length*solution.x[problem.number_variables + j];
	}
	std::cout << "New λ: "; print_vector(std::cout, this->constraint_multipliers);
	
	for (unsigned int j = 0; j < bound_multipliers.size(); j++) {
		this->bound_multipliers[j] += dual_length*delta_z[j];
	}
	/* reset z */
	current_multiplier = 0;
	for (int i = 0; i < problem.number_variables; i++) {
		if (variable_status[i] == BOUNDED_LOWER || variable_status[i] == BOUNDED_BOTH_SIDES) {
			std::max(std::min(this->bound_multipliers[current_multiplier], this->k_sigma*this->mu/trial_primal[i]), this->mu/(this->k_sigma*trial_primal[i]));
			current_multiplier++;
		}
		if (variable_status[i] == BOUNDED_UPPER || variable_status[i] == BOUNDED_BOTH_SIDES) {
			std::max(std::min(this->bound_multipliers[current_multiplier], this->k_sigma*this->mu/trial_primal[i]), this->mu/(this->k_sigma*trial_primal[i]));
			current_multiplier++;
		}
	}
	std::cout << "New z: "; print_vector(std::cout, this->bound_multipliers);
	
	/*******************************/
	/* update of barrier parameter */
	/*******************************/
	this->mu = this->update_barrier_parameter(problem, current_iterate);
	std::cout << "NEW mu VALUE: " << this->mu << "\n";
	
	LocalSolution sol(trial_primal, 0, 0);
	return sol;
}

std::vector<double> InteriorPoint::estimate_initial_multipliers(Problem& problem, Iterate& current_iterate) {
	// TODO handle number of multipliers
	
	/* build the symmetric matrix */
	COOMatrix matrix(problem.number_variables + problem.number_constraints, 0);
	
	/* identity block */
	for (int i = 0; i < problem.number_variables; i++) {
		matrix.add_term(-1., i, i);
	}
	/* Jacobian of general constraints */
	int current_column = problem.number_variables;
	for (int j = 0; j < problem.number_constraints; j++) {
		for (std::map<int,double>::iterator it = current_iterate.constraints_jacobian[j].begin(); it != current_iterate.constraints_jacobian[j].end(); it++) {
			int variable_index = it->first;
			double derivative = it->second;
			matrix.add_term(derivative, variable_index, current_column);
		}
		current_column++;
	}
	
	/* generate the right-hand side */
	std::vector<double> rhs(problem.number_variables + problem.number_constraints);
	
	/* objective gradient */
	for (std::map<int,double>::iterator it = current_iterate.objective_gradient.begin(); it != current_iterate.objective_gradient.end(); it++) {
		int variable_index = it->first;
		double derivative = it->second;
		rhs[variable_index] += problem.objective_sign*derivative;
	}	
	/* bound constraints */
	int current_multiplier = 0;
	for (int i = 0; i < problem.number_variables; i++) {
		if (this->variable_status[i] == BOUNDED_LOWER || this->variable_status[i] == BOUNDED_BOTH_SIDES) {
			rhs[i] -= this->bound_multipliers[current_multiplier];
			current_multiplier++;
		}
		if (this->variable_status[i] == BOUNDED_UPPER || this->variable_status[i] == BOUNDED_BOTH_SIDES) {
			rhs[i] -= this->bound_multipliers[current_multiplier];
			current_multiplier++;
		}
	}
	
	MA57Data data = this->solver.factorize(matrix);
	LocalSolution solution = this->solver.solve(matrix, rhs, data);
	
	/* retrieve multipliers */
	std::vector<double> multipliers(problem.number_constraints);
	for (int j = 0; j < problem.number_constraints; j++) {
		multipliers[j] = solution.x[problem.number_variables + j];
	}
	return multipliers;
}

double InteriorPoint::project_variable_in_bounds(double current_value, double lb, double ub) {
	double k1 = 1e-2;
	double k2 = 1e-2;
	
	double perturbation_lb = std::min(k1*std::max(1., std::abs(lb)), k2*(ub - lb));
	double perturbation_ub = std::min(k1*std::max(1., std::abs(ub)), k2*(ub - lb));
	current_value = std::max(current_value, lb + perturbation_lb);
	current_value = std::min(current_value, ub - perturbation_ub);
	return current_value;
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

COOMatrix InteriorPoint::generate_kkt_matrix(Problem& problem, Iterate& current_iterate, std::vector<double>& original_multipliers, std::vector<double>& variable_lb, std::vector<double>& variable_ub) {
	current_iterate.compute_hessian(problem, problem.objective_sign, original_multipliers);
	COOMatrix kkt_matrix = current_iterate.hessian.to_COO();
	kkt_matrix.size += this->constraint_multipliers.size();
	
	/* bound constraints */
	int current_multiplier = 0;
	for (int i = 0; i < problem.number_variables; i++) {
		if (this->variable_status[i] == BOUNDED_LOWER || this->variable_status[i] == BOUNDED_BOTH_SIDES) {
			kkt_matrix.add_term(this->bound_multipliers[current_multiplier]/(current_iterate.x[i] - variable_lb[i]), i, i);
			current_multiplier++;
		}
		if (this->variable_status[i] == BOUNDED_UPPER || this->variable_status[i] == BOUNDED_BOTH_SIDES) {
			kkt_matrix.add_term(this->bound_multipliers[current_multiplier]/(current_iterate.x[i] - variable_ub[i]), i, i);
			current_multiplier++;
		}
	}
	
	int current_column = problem.number_variables;
	
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
	
	/* constraint diagonal terms -L^-1 S */
	int current_slack = 0;
	current_multiplier = 0;
	for (int j = 0; j < problem.number_constraints; j++) {
		if (problem.constraint_status[j] == EQUAL_BOUNDS) {
			// no slack
			current_multiplier++;
		}
		else {
			if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
				current_column = problem.number_variables + current_multiplier;
				kkt_matrix.add_term(-(current_iterate.x[problem.number_variables + current_slack] - problem.constraint_lb[j])/this->constraint_multipliers[current_multiplier], current_column, current_column);
				current_multiplier++;
			}
			if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
				current_column = problem.number_variables + current_multiplier;
				kkt_matrix.add_term(-(current_iterate.x[problem.number_variables + current_slack] - problem.constraint_ub[j])/this->constraint_multipliers[current_multiplier], current_column, current_column);
				current_multiplier++;
			}
			current_slack++;
		}
	}
	
	/* Inertia correction */
	int current_matrix_size = kkt_matrix.matrix.size();
	for (int i = 0; i < problem.number_variables; i++) {
		kkt_matrix.add_term(this->inertia_term, i, i);
	}
	
	bool good_inertia = false;
	while (!good_inertia) {
		std::cout << "Testing factorization with inertia term " << this->inertia_term << "\n";
		this->data = this->solver.factorize(kkt_matrix);
		if (this->solver.number_negative_eigenvalues() == problem.number_constraints) {
			good_inertia = true;
			std::cout << "Factorization was a success\n";
		}
		else {
			std::cout << "Bad inertia with inertia term " << this->inertia_term << "\n";
			/* increase inertia factor */
			this->inertia_term = (this->inertia_term == 0.) ? 1e-4 : 100*this->inertia_term;
			for (int i = 0; i < problem.number_variables; i++) {
				kkt_matrix.matrix[current_matrix_size + i] = this->inertia_term;
			}
		}
	}
	
	std::cout << "KKT matrix:\n";
	for (unsigned int k = 0; k < kkt_matrix.matrix.size(); k++) {
		std::cout << "m(" << kkt_matrix.row_indices[k] << ", " << kkt_matrix.column_indices[k] << ") = " << kkt_matrix.matrix[k] << "\n";
	}
	
	return kkt_matrix;
}

std::vector<double> InteriorPoint::generate_rhs(Problem& problem, Iterate& current_iterate, std::vector<double>& original_multipliers,
		std::vector<double>& variable_lb, std::vector<double>& variable_ub) {
	/* generate the right-hand side */
	std::vector<double> rhs(problem.number_variables + this->constraint_multipliers.size());
	
	/* objective gradient */
	for (std::map<int,double>::iterator it = current_iterate.objective_gradient.begin(); it != current_iterate.objective_gradient.end(); it++) {
		int variable_index = it->first;
		double derivative = it->second;
		rhs[variable_index] -= problem.objective_sign*derivative;
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
	
	/* bound constraints */
	for (int i = 0; i < problem.number_variables; i++) {
		if (variable_status[i] == BOUNDED_LOWER || variable_status[i] == BOUNDED_BOTH_SIDES) {
			rhs[i] += this->mu/(current_iterate.x[i] - variable_lb[i]);
		}
		if (variable_status[i] == BOUNDED_UPPER || variable_status[i] == BOUNDED_BOTH_SIDES) {
			rhs[i] += this->mu/(current_iterate.x[i] - variable_ub[i]);
		}
	}
	
	int current_multiplier = 0;
	for (int j = 0; j < problem.number_constraints; j++) {
		if (problem.constraint_status[j] == EQUAL_BOUNDS) {
			rhs[problem.number_variables + current_multiplier] = -current_iterate.constraints[j];
			current_multiplier++;
		}
		if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			rhs[problem.number_variables + current_multiplier] = -current_iterate.constraints[j] + problem.constraint_lb[j] + this->mu/this->constraint_multipliers[current_multiplier];
			current_multiplier++;
		}
		if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
			rhs[problem.number_variables + current_multiplier] = -current_iterate.constraints[j] + problem.constraint_ub[j] + this->mu/this->constraint_multipliers[current_multiplier];
			current_multiplier++;
		}
	}
	std::cout << "RHS: "; print_vector(std::cout, rhs);
	
	return rhs;
}

std::vector<double> InteriorPoint::compute_slack_displacements(Problem& problem, Iterate& current_iterate, std::vector<double>& solution) {
	std::vector<double> delta_s(this->number_slacks);
	int current_slack = 0;
	for (int j = 0; j < problem.number_constraints; j++) {
		if (problem.constraint_status[j] != EQUAL_BOUNDS) {
			delta_s[current_slack] = dot(current_iterate.x, current_iterate.constraints_jacobian[j]) + current_iterate.constraints[j] - current_iterate.x[problem.number_variables + current_slack];
			current_slack++;
		}
	}
	return delta_s;
}

std::vector<double> InteriorPoint::compute_bound_multiplier_displacements(Problem& problem, Iterate& current_iterate, std::vector<double>& solution, std::vector<double>& variable_lb, std::vector<double>& variable_ub) {
	std::vector<double> delta_z(this->bound_multipliers.size());
	int current_multiplier = 0;
	for (int i = 0; i < problem.number_variables; i++) {
		if (variable_status[i] == BOUNDED_LOWER || variable_status[i] == BOUNDED_BOTH_SIDES) {
			delta_z[current_multiplier] = this->mu/(current_iterate.x[i] - variable_lb[i]) - this->bound_multipliers[current_multiplier] - this->bound_multipliers[current_multiplier]/(current_iterate.x[i] - variable_lb[i])*solution[i];
			current_multiplier++;
		}
		if (variable_status[i] == BOUNDED_UPPER || variable_status[i] == BOUNDED_BOTH_SIDES) {
			delta_z[current_multiplier] = this->mu/(current_iterate.x[i] - variable_ub[i]) - this->bound_multipliers[current_multiplier] - this->bound_multipliers[current_multiplier]/(current_iterate.x[i] - variable_ub[i])*solution[i];
			current_multiplier++;
		}
	}
	return delta_z;
}

double InteriorPoint::update_barrier_parameter(Problem& problem, Iterate& current_iterate) {
	return this->mu/10.;
}
