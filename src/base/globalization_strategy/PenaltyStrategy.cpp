#include <cmath>
#include "PenaltyStrategy.hpp"
#include "Logger.hpp"

/*
* Infeasibility detection and SQP methods for nonlinear optimization 
* http://epubs.siam.org/doi/pdf/10.1137/080738222
*/

PenaltyStrategy::PenaltyStrategy(LocalApproximation& local_approximation, double tolerance):
		GlobalizationStrategy(local_approximation, tolerance), penalty_parameter(1.) {
	this->tau = 0.5;
	this->eta = 1e-8;
	this->epsilon1 = 0.1;
	this->epsilon2 = 0.1;
}

LocalSolution PenaltyStrategy::compute_step(Problem& problem, Iterate& current_iterate, double radius) {
	/* stage a: compute the step within trust region */
	LocalSolution solution = this->local_approximation.compute_l1_penalty_step(problem, current_iterate, radius, this->penalty_parameter, this->penalty_constraints);
	DEBUG << solution;
	
	/* if penalty parameter is already 0, no need to decrease it */
	if (0. < this->penalty_parameter) {
		/* check infeasibility */
		double linear_model = this->compute_linear_model(problem, solution);
		if (linear_model != 0.) {
			double current_penalty_parameter = this->penalty_parameter;
			
			/* stage c: solve the ideal l1 penalty problem with a zero penalty (no objective) */
			LocalSolution ideal_solution = this->local_approximation.compute_l1_penalty_step(problem, current_iterate, radius, 0., this->penalty_constraints);
			DEBUG << ideal_solution;
			
			/* stage f: update the penalty parameter */
			std::vector<double> ideal_multipliers = this->compute_multipliers(problem, ideal_solution);
			/* compute the ideal error (with a zero penalty parameter) */
			double ideal_error = this->compute_error(problem, current_iterate, ideal_multipliers, 0.);
			
			if (ideal_error == 0.) {
				/* stage f: update the penalty parameter */
				this->penalty_parameter = 0.;
			}
			else {
				double ideal_linear_model = this->compute_linear_model(problem, ideal_solution);
				
				/* decrease penalty parameter to satisfy 2 conditions */
				bool condition1 = false, condition2 = false;
				while (!condition2) {
					this->penalty_parameter *= this->tau;
					if (this->penalty_parameter < 1e-10) {
						this->penalty_parameter = 0.;
						condition2 = true;
					}
					
					DEBUG << "\nSolving with penalty parameter " << this->penalty_parameter << "\n";
					solution = this->local_approximation.compute_l1_penalty_step(problem, current_iterate, radius, this->penalty_parameter, this->penalty_constraints);
					DEBUG << solution;
					
					double trial_linear_model = this->compute_linear_model(problem, solution);
					if (!condition1) {
						/* stage d: reach a fraction of the ideal decrease */
						if((ideal_linear_model == 0. && trial_linear_model == 0.) || (ideal_linear_model != 0. &&
							current_iterate.residual - trial_linear_model >= this->epsilon1*(current_iterate.residual - ideal_linear_model))) {
							condition1 = true;
						}
					}
					/* stage e: further decrease penalty parameter if necessary */
					if (condition1 && current_iterate.residual - solution.objective >= this->epsilon2*(current_iterate.residual - ideal_solution.objective)) {
						condition2 = true;
					}
				}
				
				/* stage f: update the penalty parameter */
				double term = ideal_error / std::max(1., current_iterate.residual);
				this->penalty_parameter = std::min(this->penalty_parameter, term*term);
			}
			
			if (this->penalty_parameter < current_penalty_parameter) {
				DEBUG << "Penalty parameter updated to " << this->penalty_parameter << "\n";
				/* recompute the solution */
				if (this->penalty_parameter == 0.) {
					solution = ideal_solution;
				}
				else {
					solution = this->local_approximation.compute_l1_penalty_step(problem, current_iterate, radius, this->penalty_parameter, this->penalty_constraints);
					DEBUG << solution;
				}
			}
		}
	}
	INFO << "penalty parameter: " << this->penalty_parameter << "\t";
	return solution;
}

bool PenaltyStrategy::check_step(Problem& problem, Iterate& current_iterate, LocalSolution& solution, double step_length) {
	/* stage g: line-search along fixed step */
	
	/* retrieve only original primal and dual variables from the step */
	std::vector<double> d(problem.number_variables);
	for (int i = 0; i < problem.number_variables; i++) {
		d[i] = solution.x[i];
	}
	std::vector<double> multipliers = this->compute_multipliers(problem, solution);
	
	/* generate the trial point */
	std::vector<double> x_trial = add_vectors(current_iterate.x, d, step_length);
	Iterate trial_iterate(problem, x_trial, multipliers);
	
	/* compute current exact l1 penalty: rho f + sum max(0, c) */
	double current_exact_l1_penalty = this->penalty_parameter*current_iterate.objective + current_iterate.residual;
	/* compute trial exact l1 penalty */
	double trial_exact_l1_penalty = this->penalty_parameter*trial_iterate.objective + trial_iterate.residual;
	
	/* check the validity of the trial step */
	bool accept = false;
	if (current_exact_l1_penalty - trial_exact_l1_penalty >= this->eta*step_length*(current_iterate.residual - solution.objective)) {
		accept = true;
		trial_iterate.KKTerror = this->compute_KKT_error(problem, trial_iterate);
		trial_iterate.complementarity_error = this->compute_complementarity_error(problem, trial_iterate);
		double step_norm = step_length*norm_inf(d);
		trial_iterate.status = this->compute_status(problem, trial_iterate, step_norm);
		current_iterate = trial_iterate;
	}
	return accept;
}

OptimalityStatus PenaltyStrategy::compute_status(Problem& problem, Iterate& trial_iterate, double step_norm) {
	OptimalityStatus status = NOT_OPTIMAL;
	
	/* test for optimality */
	double optimality_error = this->compute_error(problem, trial_iterate, trial_iterate.multipliers, this->penalty_parameter);
	DEBUG << "Ek(lambda_k, rho_k) = " << optimality_error << "\n";
	if (optimality_error <= this->tolerance && trial_iterate.residual <= this->tolerance*problem.number_constraints) {
		status = KKT_POINT;
		/* rescale the multipliers */
		if (0. < this->penalty_parameter) {
			for (unsigned int k = 0; k < trial_iterate.multipliers.size(); k++) {
				trial_iterate.multipliers[k] /= this->penalty_parameter;
			}
		}
	}
	else {
		double infeasibility_error = this->compute_error(problem, trial_iterate, trial_iterate.multipliers, 0.);
		DEBUG << "Ek(lambda_k, 0.) = " << infeasibility_error << "\n";
		if (infeasibility_error <= this->tolerance && trial_iterate.residual > this->tolerance*problem.number_constraints) {
			status = FJ_POINT;
		}
		else if (step_norm <= this->tolerance/100.) {
			if (trial_iterate.residual <= this->tolerance*problem.number_constraints) {
				status = FEASIBLE_SMALL_STEP;
			}
			else {
				status = INFEASIBLE_SMALL_STEP;
			}
		}
	}
	return status;
}

double PenaltyStrategy::compute_linear_model(Problem& problem, LocalSolution& solution) {
	double linear_model = 0.;
	for (int k = 0; k < this->penalty_constraints.number_additional_variables; k++) {
		linear_model += solution.x[problem.number_variables + k];
	}
	return linear_model;
}

std::vector<double> PenaltyStrategy::compute_multipliers(Problem& problem, LocalSolution& solution) {
	std::vector<double> multipliers(problem.number_variables + problem.number_constraints);
	for (int i = 0; i < problem.number_variables; i++) {
		multipliers[i] = solution.multipliers[i];
	}
	int current_constraint = problem.number_variables + this->penalty_constraints.number_additional_variables;
	for (int j = 0; j < problem.number_constraints; j++) {
		if (this->penalty_constraints.status[j] == BOUNDED_BOTH_SIDES) {
			/* only case where 2 constraints were generated */
			/* only one bound is active: one multiplier is > 0, the other is 0 */
			multipliers[problem.number_variables + j] = solution.multipliers[current_constraint] + solution.multipliers[current_constraint+1];
			current_constraint += 2;
		}
		else {
			/* only 1 constraint was generated */
			multipliers[problem.number_variables + j] = solution.multipliers[current_constraint];
			current_constraint++;
		}
	}
	return multipliers;
}

double PenaltyStrategy::compute_error(Problem& problem, Iterate& current_iterate, std::vector<double>& multipliers, double penalty_parameter) {
	/* measure that combines KKT error and complementarity error */
	double error = 0.;
	
	/* KKT error */
	std::vector<double> lagrangian_gradient = this->compute_lagrangian_gradient(problem, current_iterate, penalty_parameter, multipliers);
	/* compute 1-norm */
	error += norm_1(lagrangian_gradient);
		
	/* complementarity error */
	
	/* bound constraints */
	for (int i = 0; i < problem.number_variables; i++) {
		if (problem.variable_lb[i] < current_iterate.x[i] && current_iterate.x[i] < problem.variable_ub[i]) {
			double multiplier_i = multipliers[i];
			
			if (multiplier_i > 0.) {
				error += std::abs(multiplier_i*(current_iterate.x[i] - problem.variable_lb[i]));
			}
			else if (multiplier_i < 0.) {
				error += std::abs(multiplier_i*(current_iterate.x[i] - problem.variable_ub[i]));
			}
		}
	}
	/* check if constraint is strictly satisfied or violated */
	for (int j = 0; j < problem.number_constraints; j++) {
		double multiplier_j = multipliers[problem.number_variables + j];
		
		/* violated */
		if (current_iterate.constraints[j] < problem.constraint_lb[j]) {
			error += std::abs((1. - multiplier_j)*(current_iterate.constraints[j] - problem.constraint_lb[j]));
		}
		else if (problem.constraint_ub[j] < current_iterate.constraints[j]) {
			error += std::abs((1. - multiplier_j)*(current_iterate.constraints[j] - problem.constraint_ub[j]));
			
		}
		else {
			/* active or strictly satisfied */
			if (multiplier_j > 0.) {
				error += std::abs(multiplier_j*(current_iterate.constraints[j] - problem.constraint_lb[j]));
			}
			else if (multiplier_j < 0.) {
				error += std::abs(multiplier_j*(current_iterate.constraints[j] - problem.constraint_ub[j]));
			}
		}
	}
	return error;
}

	//for (int j = 0; j < problem.number_constraints; j++) {
		//double multiplier_j = multipliers[problem.number_variables + j];
		///* strictly satisfied */
		//if (problem.constraint_lb[j] < current_iterate.constraints[j] && current_iterate.constraints[j] < problem.constraint_ub[j]) {
			//if (multiplier_j > 0.) {
				//error += std::abs(multiplier_j*(current_iterate.constraints[j] - problem.constraint_lb[j]));
			//}
			//else if (multiplier_j < 0.) {
				//error += std::abs(multiplier_j*(current_iterate.constraints[j] - problem.constraint_ub[j]));
			//}
		//}
		///* violated */
		//else if (current_iterate.constraints[j] < problem.constraint_lb[j]) {
			//error += std::abs((1. - multiplier_j)*(current_iterate.constraints[j] - problem.constraint_lb[j]));
		//}
		//else if (problem.constraint_ub[j] < current_iterate.constraints[j]) {
			//error += std::abs((1. - multiplier_j)*(current_iterate.constraints[j] - problem.constraint_ub[j]));
			
		//}
	//}

double PenaltyStrategy::compute_KKT_error(Problem& problem, Iterate& current_iterate) {
	std::vector<double> lagrangian_gradient = this->compute_lagrangian_gradient(problem, current_iterate, this->penalty_parameter, current_iterate.multipliers);
	double KKTerror = norm_2(lagrangian_gradient);
	return KKTerror;
}

void PenaltyStrategy::initialize(Problem& problem, Iterate& current_iterate) {
	/* compute the number of necessary additional variables and constraints */
	this->penalty_constraints.status.resize(problem.number_constraints);
	this->penalty_constraints.number_additional_variables = 0;
	this->penalty_constraints.number_constraints = 0;
	
	for (int j = 0; j < problem.number_constraints; j++) {
		if (problem.constraint_lb[j] == problem.constraint_ub[j]) {
			this->penalty_constraints.status[j] = EQUALITY;
			this->penalty_constraints.number_additional_variables += 2;
			this->penalty_constraints.number_constraints++;
		}
		else if (-INFINITY < problem.constraint_lb[j] && problem.constraint_ub[j] < INFINITY) {
			this->penalty_constraints.status[j] = BOUNDED_BOTH_SIDES;
			this->penalty_constraints.number_additional_variables += 2;
			this->penalty_constraints.number_constraints += 2;
		}
		else if (-INFINITY < problem.constraint_lb[j]) {
			this->penalty_constraints.status[j] = BOUNDED_LOWER;
			this->penalty_constraints.number_additional_variables++;
			this->penalty_constraints.number_constraints++;
		}
		else {
			this->penalty_constraints.status[j] = BOUNDED_UPPER;
			this->penalty_constraints.number_additional_variables++;
			this->penalty_constraints.number_constraints++;
		}
	}
	
	/* allocate the subproblem solver */
	int number_variables = problem.number_variables + this->penalty_constraints.number_additional_variables;
	int number_constraints = this->penalty_constraints.number_constraints;
	this->local_approximation.allocate_solver(number_variables, number_constraints);
	return;
}
