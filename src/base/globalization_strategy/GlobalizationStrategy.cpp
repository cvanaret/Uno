#include <ostream>
#include <cmath>
#include "GlobalizationStrategy.hpp"

GlobalizationStrategy::GlobalizationStrategy(Subproblem& subproblem, double tolerance) : subproblem(subproblem) {
    this->tolerance = tolerance;
}

GlobalizationStrategy::~GlobalizationStrategy() {
}

std::vector<double> GlobalizationStrategy::compute_lagrangian_gradient(Problem& problem, Iterate& current_iterate, double objective_multiplier, Multipliers& multipliers) {
    std::vector<double> lagrangian_gradient(problem.number_variables);

    /* objective gradient */
    if (objective_multiplier != 0.) {
        current_iterate.compute_objective_gradient(problem);
        
        for (std::pair<int, double> term : current_iterate.objective_gradient) {
            int variable_index = term.first;
            double derivative = term.second;
            /* scale the objective gradient */
            lagrangian_gradient[variable_index] += objective_multiplier*derivative;
        }
    }
    /* bound constraints */
    for (int i = 0; i < problem.number_variables; i++) {
        lagrangian_gradient[i] -= multipliers.bounds[i];
    }
    /* constraints */
    current_iterate.compute_constraints_jacobian(problem);

    for (int j = 0; j < problem.number_constraints; j++) {
        double multiplier_j = multipliers.constraints[j];
        if (multiplier_j != 0.) {
            for (std::pair<int, double> term : current_iterate.constraints_jacobian[j]) {
                int variable_index = term.first;
                double derivative = term.second;
                lagrangian_gradient[variable_index] -= multiplier_j*derivative;
            }
        }
    }
    return lagrangian_gradient;
}
/* complementary slackness error. Use abs/1e-8 to safeguard */
double GlobalizationStrategy::compute_complementarity_error(const Problem& problem, Iterate& iterate) {
    double complementarity_error = 0.;
	
    /* bound constraints */
    for (int i = 0; i < problem.number_variables; i++) {
        double multiplier_i = iterate.multipliers.bounds[i];

        if (multiplier_i > this->tolerance / 10.) {
            complementarity_error += std::abs(multiplier_i * (iterate.x[i] - problem.variables_bounds[i].lb));
        }
        else if (multiplier_i < -this->tolerance / 10.) {
            complementarity_error += std::abs(multiplier_i * (iterate.x[i] - problem.variables_bounds[i].ub));
        }
    }
    /* constraints */
    for (int j = 0; j < problem.number_constraints; j++) {
        double multiplier_j = iterate.multipliers.constraints[j];

        if (multiplier_j > this->tolerance / 10.) {
            complementarity_error += std::abs(multiplier_j * (iterate.constraints[j] - problem.constraints_bounds[j].lb));
        }
        else if (multiplier_j < -this->tolerance / 10.) {
            complementarity_error += std::abs(multiplier_j * (iterate.constraints[j] - problem.constraints_bounds[j].ub));
        }
    }
    return complementarity_error;
}

double GlobalizationStrategy::compute_KKT_error(Problem& problem, Iterate& iterate, double objective_mutiplier) {
    std::vector<double> lagrangian_gradient = this->compute_lagrangian_gradient(problem, iterate, objective_mutiplier, iterate.multipliers);
    double KKTerror = norm_2(lagrangian_gradient);
    return KKTerror;
}