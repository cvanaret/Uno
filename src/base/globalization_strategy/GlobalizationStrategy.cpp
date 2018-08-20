#include <ostream>
#include <cmath>
#include "GlobalizationStrategy.hpp"

GlobalizationStrategy::GlobalizationStrategy(Subproblem& subproblem, double tolerance) : subproblem(subproblem) {
    this->tolerance = tolerance;
}

GlobalizationStrategy::~GlobalizationStrategy() {
}

std::vector<double> GlobalizationStrategy::compute_lagrangian_gradient(Problem& problem, Iterate& current_iterate, double objective_multiplier,
                                                                       std::vector<double>& bound_multipliers, std::vector<double>& constraint_multipliers) {
    std::vector<double> lagrangian_gradient(problem.number_variables);

    /* objective gradient */
    if (objective_multiplier != 0.) {
        if (!current_iterate.is_objective_gradient_computed) {
            std::map<int, double> objective_gradient = problem.objective_sparse_gradient(current_iterate.x);
            current_iterate.set_objective_gradient(objective_gradient);
        }
        for (std::pair<int, double> term : current_iterate.objective_gradient) {
            int variable_index = term.first;
            double derivative = term.second;
            /* scale the objective gradient */
            lagrangian_gradient[variable_index] += objective_multiplier*derivative;
        }
    }
    /* bound constraints */
    for (int i = 0; i < problem.number_variables; i++) {
        lagrangian_gradient[i] -= bound_multipliers[i];
    }
    /* constraints */
    current_iterate.compute_constraints_jacobian(problem);

    for (int j = 0; j < problem.number_constraints; j++) {
        double multiplier_j = constraint_multipliers[j];
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
double GlobalizationStrategy::compute_complementarity_error(const Problem& problem, Iterate& current_iterate) {
    double complementarity_error = 0.;

    std::cout << "complementarity error = 0";

    /* bound constraints */
    for (int i = 0; i < problem.number_variables; i++) {
        double multiplier_i = current_iterate.bound_multipliers[i];

        if (multiplier_i > this->tolerance / 10.) {
            std::cout << " + abs(" << multiplier_i << "*(x[" << i << "] - " << problem.variable_lb[i] << "))";
            complementarity_error += std::abs(multiplier_i * (current_iterate.x[i] - problem.variable_lb[i]));
        }
        else if (multiplier_i < -this->tolerance / 10.) {
            std::cout << " + abs(" << multiplier_i << "*(x[" << i << "] - " << problem.variable_ub[i] << "))";
            complementarity_error += std::abs(multiplier_i * (current_iterate.x[i] - problem.variable_ub[i]));
        }
    }
    /* constraints */
    for (int j = 0; j < problem.number_constraints; j++) {
        double multiplier_j = current_iterate.constraint_multipliers[j];

        if (multiplier_j > this->tolerance / 10.) {
            complementarity_error += std::abs(multiplier_j * (current_iterate.constraints[j] - problem.constraint_lb[j]));
        }
        else if (multiplier_j < -this->tolerance / 10.) {
            complementarity_error += std::abs(multiplier_j * (current_iterate.constraints[j] - problem.constraint_ub[j]));
        }
    }
    std::cout << "\n";
    return complementarity_error;
}
