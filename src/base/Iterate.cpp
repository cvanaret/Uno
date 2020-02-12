#include "Iterate.hpp"
#include "Utils.hpp"
#include "Logger.hpp"

Iterate::Iterate(Problem& problem, std::vector<double>& x, Multipliers& multipliers) :
x(x), multipliers(multipliers) {

    /* project point onto bounds */
    for (int i = 0; i < problem.number_variables; i++) {
        this->x[i] = std::max(this->x[i], problem.variable_lb[i]);
        this->x[i] = std::min(this->x[i], problem.variable_ub[i]);
    }

    /* objective */
    this->objective = problem.objective(this->x);

    /* constraints */
    this->constraints = problem.evaluate_constraints(this->x);
    this->residual = problem.l1_norm(this->constraints);
    //this->feasibility_measure = this->residual;
    //this->optimality_measure = this->objective;

    /* Jacobian and Hessian will be evaluated only when necessary */
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

std::ostream& operator<<(std::ostream &stream, Iterate& iterate) {
    stream << "x: ";
    print_vector(stream, iterate.x, 0, 50);

    stream << "Bound multipliers: ";
    print_vector(stream, iterate.multipliers.bounds, 0, 50);

    stream << "Constraint multipliers: ";
    print_vector(stream, iterate.multipliers.constraints, 0, 50);

    stream << "Objective: " << iterate.objective << "\n";

    //stream << "Constraints:";
    //for (double cj: iterate.constraints) {
    //	stream << " " << cj;
    //}
    //stream << "\n";

    stream << "Residual: " << iterate.residual << "\n";

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
