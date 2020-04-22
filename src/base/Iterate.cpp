#include "Iterate.hpp"
#include "Utils.hpp"
#include "Logger.hpp"

Iterate::Iterate(Problem& problem, std::vector<double>& x, Multipliers& multipliers, int residual_norm) : x(x), multipliers(multipliers) {
    /* objective */
    // TODO: do not evaluate
    this->objective = problem.objective(this->x);

    /* constraints */
    this->constraints = problem.evaluate_constraints(this->x);
    this->residual = problem.infeasible_residual_norm(this->constraints, residual_norm);
    //this->feasibility_measure = this->residual;
    //this->optimality_measure = this->objective;

    /* evaluations will be performed only when necessary */
    this->is_objective_gradient_computed = false;
    this->is_constraints_jacobian_computed = false;
    this->is_hessian_computed = false;

    /* status */
    this->status = NOT_OPTIMAL;
}

void Iterate::set_objective_gradient(std::map<int, double>& objective_gradient) {
    this->objective_gradient = objective_gradient;
    this->is_objective_gradient_computed = true;
    return;
}

void Iterate::compute_objective_gradient(Problem& problem) {
    if (!this->is_objective_gradient_computed) {
        this->objective_gradient = problem.objective_sparse_gradient(this->x);
        this->is_objective_gradient_computed = true;
    }
    return;
}

void Iterate::compute_constraints_jacobian(Problem& problem) {
    if (!this->is_constraints_jacobian_computed) {
        this->constraints_jacobian = problem.constraints_sparse_jacobian(this->x);
        this->is_constraints_jacobian_computed = true;
    }
    return;
}

void Iterate::compute_hessian(Problem& problem, double objective_multiplier, std::vector<double>& constraint_multipliers) {
    if (!this->is_hessian_computed) {
        this->hessian = problem.lagrangian_hessian(this->x, objective_multiplier, constraint_multipliers);
        this->is_hessian_computed = true;
    }
    return;
}

std::vector<double> Iterate::lagrangian_gradient(Problem& problem) {
    std::vector<double> lagrangian_gradient(this->x.size());

    /* objective gradient */
    if (problem.objective_sign != 0.) {
        this->compute_objective_gradient(problem);
        
        /* scale the objective gradient */
        for (std::pair<int, double> term : this->objective_gradient) {
            int variable_index = term.first;
            double derivative = term.second;
            lagrangian_gradient[variable_index] += problem.objective_sign*derivative;
        }
    }
    /* bound constraints */
    for (unsigned int i = 0; i < this->x.size(); i++) {
        lagrangian_gradient[i] += - this->multipliers.lower_bounds[i] - this->multipliers.upper_bounds[i];
    }
    /* constraints */
    this->compute_constraints_jacobian(problem);

    for (int j = 0; j < problem.number_constraints; j++) {
        double multiplier_j = this->multipliers.constraints[j];
        if (multiplier_j != 0.) {
            for (std::pair<int, double> term : this->constraints_jacobian[j]) {
                int variable_index = term.first;
                double derivative = term.second;
                lagrangian_gradient[variable_index] -= multiplier_j*derivative;
            }
        }
    }
    return lagrangian_gradient;
}

std::ostream& operator<<(std::ostream &stream, Iterate& iterate) {
    stream << "x: "; print_vector(stream, iterate.x, 0, 50);
    stream << "Lower bound multipliers: "; print_vector(stream, iterate.multipliers.lower_bounds, 0, 50);
    stream << "Upper bound multipliers: "; print_vector(stream, iterate.multipliers.upper_bounds, 0, 50);
    stream << "Constraint multipliers: "; print_vector(stream, iterate.multipliers.constraints, 0, 50);
    stream << "Objective value: " << iterate.objective << "\n";

    //stream << "Constraints:";
    //for (double cj: iterate.constraints) {
    //	stream << " " << cj;
    //}
    //stream << "\n";

    stream << "Constraint residual: " << iterate.residual << "\n";

    stream << "Optimality measure: " << iterate.optimality_measure << "\n";
    stream << "Feasibility measure: " << iterate.feasibility_measure << "\n";

    stream << "KKT error: " << iterate.KKTerror << "\n";
    stream << "Complementarity error: " << iterate.complementarity_error << "\n";
    return stream;
}

std::ostream& operator<<(std::ostream &stream, OptimalityStatus& status) {
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
