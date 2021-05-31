#include <limits>
#include "Iterate.hpp"
#include "Vector.hpp"
#include "Logger.hpp"

int Iterate::number_eval_objective = 0;
int Iterate::number_eval_constraints = 0;
int Iterate::number_eval_jacobian = 0;
int Iterate::number_eval_hessian = 0;

Iterate::Iterate(const std::vector<double>& x, const Multipliers& multipliers) : x(x), multipliers(multipliers),
    objective(std::numeric_limits<double>::infinity()), is_objective_computed(false), are_constraints_computed(false), is_objective_gradient_computed(false), is_constraints_jacobian_computed(false), hessian(x.size(), 1), is_hessian_computed(false),
    status(NOT_OPTIMAL),
    residuals({0., 0., 0., 0.}), feasibility_measure(0.), optimality_measure(0.) {
}

//std::optional<Iterate>{trial_iterate};
//return std::nullopt;

void Iterate::compute_objective(const Problem& problem) {
    if (!this->is_objective_computed) {
        this->objective = problem.objective(this->x);
        this->is_objective_computed = true;
        Iterate::number_eval_objective++;
    }
    return;
}

void Iterate::compute_constraints(const Problem& problem) {
    if (!this->are_constraints_computed) {
        this->constraints = problem.evaluate_constraints(this->x);
        this->are_constraints_computed = true;
        Iterate::number_eval_constraints++;
    }
    return;
}

void Iterate::compute_objective_gradient(const Problem& problem) {
    if (!this->is_objective_gradient_computed) {
        this->objective_gradient = problem.objective_gradient(this->x);
        this->is_objective_gradient_computed = true;
    }
    return;
}

void Iterate::set_objective_gradient(const SparseGradient& objective_gradient) {
    this->objective_gradient = objective_gradient;
    this->is_objective_gradient_computed = true;
    return;
}

void Iterate::compute_constraints_jacobian(const Problem& problem) {
    if (!this->is_constraints_jacobian_computed) {
        this->constraints_jacobian = problem.constraints_jacobian(this->x);
        this->is_constraints_jacobian_computed = true;
        Iterate::number_eval_jacobian++;
    }
    return;
}

void Iterate::compute_hessian(const Problem& problem, double objective_multiplier, const std::vector<double>& constraint_multipliers) {
    if (!this->is_hessian_computed) {
        this->hessian = problem.lagrangian_hessian(this->x, objective_multiplier, constraint_multipliers);
        this->is_hessian_computed = true;
        Iterate::number_eval_hessian++;
    }
    return;
}

std::vector<double> Iterate::lagrangian_gradient(const Problem& problem, double objective_mutiplier, const Multipliers& multipliers) {
    std::vector<double> lagrangian_gradient(problem.number_variables);

    /* objective gradient */
    if (objective_mutiplier != 0.) {
        this->compute_objective_gradient(problem);
        
        /* scale the objective gradient */
        for (std::pair<int, double> term : this->objective_gradient) {
            int i = term.first;
            double derivative = term.second;
            if (i < problem.number_variables) {
                lagrangian_gradient[i] += objective_mutiplier*derivative;
            }
        }
    }
    /* bound constraints */
    for (int i = 0; i < problem.number_variables; i++) {
        lagrangian_gradient[i] -= multipliers.lower_bounds[i] + multipliers.upper_bounds[i];
    }
    
    /* constraints */
    this->compute_constraints_jacobian(problem);
    for (int j = 0; j < problem.number_constraints; j++) {
        double multiplier_j = multipliers.constraints[j];
        if (multiplier_j != 0.) {
            for (std::pair<int, double> term : this->constraints_jacobian[j]) {
                int i = term.first;
                double derivative = term.second;
                if (i < problem.number_variables) {
                    lagrangian_gradient[i] -= multiplier_j*derivative;
                }
            }
        }
    }
    return lagrangian_gradient;
}

std::ostream& operator<<(std::ostream &stream, const Iterate& iterate) {
    stream << "x: "; print_vector(stream, iterate.x);
    stream << "Lower bound multipliers: "; print_vector(stream, iterate.multipliers.lower_bounds);
    stream << "Upper bound multipliers: "; print_vector(stream, iterate.multipliers.upper_bounds);
    stream << "Constraint multipliers: "; print_vector(stream, iterate.multipliers.constraints);
    stream << "Objective value: " << iterate.objective << "\n";

    stream << "Constraint residual: " << iterate.residuals.constraints << "\n";
    stream << "KKT residual: " << iterate.residuals.KKT << "\n";
    stream << "FJ residual: " << iterate.residuals.KKT << "\n";
    stream << "Complementarity residual: " << iterate.residuals.complementarity << "\n";
    
    stream << "Optimality measure: " << iterate.optimality_measure << "\n";
    stream << "Feasibility measure: " << iterate.feasibility_measure << "\n";
    return stream;
}

std::ostream& operator<<(std::ostream &stream, TerminationStatus status) {
    if (status == NOT_OPTIMAL) {
        stream << "not optimal";
    }
    else if (status == KKT_POINT) {
        stream << "KKT point";
    }
    else if (status == FJ_POINT) {
        stream << "FJ point";
    }
    else if (status == FEASIBLE_SMALL_STEP) {
        stream << "feasible small step";
    }
    else if (status == INFEASIBLE_SMALL_STEP) {
        stream << "infeasible small step";
    }
    return stream;
}
