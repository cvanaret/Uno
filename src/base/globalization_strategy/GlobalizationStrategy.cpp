#include <ostream>
#include <cmath>
#include "GlobalizationStrategy.hpp"

GlobalizationStrategy::GlobalizationStrategy(LocalApproximation& local_approximation, double tolerance): local_approximation(local_approximation){
	this->tolerance = tolerance;
}

std::vector<double> GlobalizationStrategy::compute_lagrangian_gradient(Problem& problem, Iterate& current_point, double objective_multiplier, std::vector<double>& multipliers) {
	std::vector<double> lagrangian_gradient(problem.number_variables);
	
	/* objective gradient */
	if (objective_multiplier != 0.) {
		if (!current_point.is_objective_gradient_computed) {
			std::map<int,double> objective_gradient = problem.objective_sparse_gradient(current_point.x);
			current_point.set_objective_gradient(objective_gradient);
		}
		for (std::map<int,double>::iterator it = current_point.objective_gradient.begin(); it != current_point.objective_gradient.end(); it++) {
			int variable_index = it->first;
			double derivative = it->second;
			/* scale the objective gradient */
			lagrangian_gradient[variable_index] += objective_multiplier*derivative;
		}
	}
	/* bound constraints */
	for (int i = 0; i < problem.number_variables; i++) {
		double multiplier_i = multipliers[i];
		lagrangian_gradient[i] -= multiplier_i;
	}
	/* constraints */
	current_point.compute_constraint_jacobian(problem);
	
	for (int j = 0; j < problem.number_constraints; j++) {
		double multiplier_j = multipliers[problem.number_variables + j];
		if (multiplier_j != 0.) {
			for (std::map<int,double>::iterator it = current_point.constraint_jacobian[j].begin(); it != current_point.constraint_jacobian[j].end(); it++) {
				int variable_index = it->first;
				double derivative = it->second;
				lagrangian_gradient[variable_index] -= multiplier_j*derivative;
			}
		}
	}
	return lagrangian_gradient;
}

/* complementary slackness error. Use abs/1e-8 to safeguard */
double GlobalizationStrategy::compute_complementarity_error(const Problem& problem, Iterate& current_point) {
	double complementarity_error = 0.;
	
	/* bound constraints */
	for (int i = 0; i < problem.number_variables; i++) {
		double multiplier_i = current_point.multipliers[i];
		
		if (multiplier_i > this->tolerance/10.) {
			complementarity_error += std::abs(multiplier_i*(current_point.x[i] - problem.variable_lb[i]));
		}
		else if (multiplier_i < -this->tolerance/10.) {
			complementarity_error += std::abs(multiplier_i*(current_point.x[i] - problem.variable_ub[i]));
		}
		
	}
	/* constraints */
	for(int j = 0; j < problem.number_constraints; j++) {
		double multiplier_j = current_point.multipliers[problem.number_variables + j];
		
		if (multiplier_j > this->tolerance/10.) {
			complementarity_error += std::abs(multiplier_j*(current_point.constraints[j] - problem.constraint_lb[j]));
		}
		else if (multiplier_j < -this->tolerance/10.) {
			complementarity_error += std::abs(multiplier_j*(current_point.constraints[j] - problem.constraint_ub[j]));
		}
	}
	return complementarity_error;
}
