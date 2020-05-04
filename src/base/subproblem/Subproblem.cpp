#include "Subproblem.hpp"

Subproblem::Subproblem(std::string restoration_norm) : residual_norm(restoration_norm), number_subproblems_solved(0), subproblem_definition_changed(false) {
}

Subproblem::~Subproblem() {
}

std::vector<Range> Subproblem::generate_variables_bounds(Iterate& current_iterate, double trust_region_radius) {
    std::vector<Range> variables_bounds(current_iterate.x.size());
    /* bounds intersected with trust region  */
    for (unsigned int i = 0; i < current_iterate.x.size(); i++) {
        double lb = std::max(-trust_region_radius, this->subproblem_variables_bounds[i].lb - current_iterate.x[i]);
        double ub = std::min(trust_region_radius, this->subproblem_variables_bounds[i].ub - current_iterate.x[i]);
        variables_bounds[i] = {lb, ub};
    }
    return variables_bounds;
}

double Subproblem::project_variable_in_bounds(double variable_value, Range& variable_bounds) {
    double k1 = 1e-2;
    double k2 = 1e-2;

    double perturbation_lb = std::min(k1 * std::max(1., std::abs(variable_bounds.lb)), k2 * (variable_bounds.ub - variable_bounds.lb));
    double perturbation_ub = std::min(k1 * std::max(1., std::abs(variable_bounds.ub)), k2 * (variable_bounds.ub - variable_bounds.lb));
    variable_value = std::max(variable_value, variable_bounds.lb + perturbation_lb);
    variable_value = std::min(variable_value, variable_bounds.ub - perturbation_ub);
    return variable_value;
}

std::vector<Range> Subproblem::generate_constraints_bounds(Problem& problem, std::vector<double>& current_constraints) {
    std::vector<Range> constraints_bounds(problem.number_constraints);
    for (int j = 0; j < problem.number_constraints; j++) {
        double lb = problem.constraints_bounds[j].lb - current_constraints[j];
        double ub = problem.constraints_bounds[j].ub - current_constraints[j];
        constraints_bounds[j] = {lb, ub};
    }
    return constraints_bounds;
}

std::vector<double> Subproblem::compute_least_square_multipliers(Problem& problem, Iterate& current_iterate, std::vector<double>& default_multipliers, double multipliers_max_size) {
    MA57Solver solver;
    return Subproblem::compute_least_square_multipliers(problem, current_iterate, default_multipliers, solver, multipliers_max_size);
}

std::vector<double> Subproblem::compute_least_square_multipliers(Problem& problem, Iterate& current_iterate, std::vector<double>& default_multipliers, MA57Solver& solver, double multipliers_max_size) {
    current_iterate.compute_objective_gradient(problem);
    current_iterate.compute_constraints_jacobian(problem);
    
    /******************************/
    /* build the symmetric matrix */
    /******************************/
    COOMatrix matrix(current_iterate.x.size() + problem.number_constraints);
    
    /* identity blocks */
    for (unsigned int i = 0; i < current_iterate.x.size(); i++) {
        matrix.add_term(1., i, i);
    }
    /* Jacobian of general constraints */
    for (int j = 0; j < problem.number_constraints; j++) {
        for (std::pair<const int, double>& term : current_iterate.constraints_jacobian[j]) {
            int variable_index = term.first;
            double derivative = term.second;
            matrix.add_term(derivative, variable_index, current_iterate.x.size() + j);
        }
    }
    DEBUG << "Multipliers estimation: KKT matrix:\n";
    for (unsigned int k = 0; k < matrix.matrix.size(); k++) {
        DEBUG << "m(" << matrix.row_indices[k] << ", " << matrix.column_indices[k] << ") = " << matrix.matrix[k] << "\n";
    }
    
    /********************************/
    /* generate the right-hand side */
    /********************************/
    std::vector<double> rhs(current_iterate.x.size() + problem.number_constraints);

    /* objective gradient */
    for (std::pair<int, double> term : current_iterate.objective_gradient) {
        int i = term.first;
        double derivative = term.second;
        rhs[i] += problem.objective_sign*derivative;
    }
    /* variable bound constraints */
    for (unsigned int i = 0; i < current_iterate.x.size(); i++) {
        rhs[i] -= current_iterate.multipliers.lower_bounds[i];
        rhs[i] -= current_iterate.multipliers.upper_bounds[i];
    }    
    DEBUG << "Multipliers RHS:\n"; print_vector(DEBUG, rhs);
    
    MA57Factorization factorization = solver.factorize(matrix);
    std::vector<double> solution = solver.solve(factorization, rhs);
    DEBUG << "Solution: "; print_vector(DEBUG, solution);
    
    /* retrieve multipliers */
    std::vector<double> multipliers(problem.number_constraints);
    for (int j = 0; j < problem.number_constraints; j++) {
        multipliers[j] = solution[current_iterate.x.size() + j];
    }
    // if multipliers too big, discard them
    if (norm_inf(multipliers) > multipliers_max_size) {
        return default_multipliers;
    }
    return multipliers;
}

/* complementary slackness error. Use abs/1e-8 to safeguard */
double Subproblem::compute_complementarity_error(Problem& problem, Iterate& iterate, Multipliers& multipliers) {
    std::cout << "Calling Subproblem::compute_complementarity_error\n";
    double complementarity_error = 0.;
    /* bound constraints */
    for (int i = 0; i < problem.number_variables; i++) {
        if (-INFINITY < problem.variables_bounds[i].lb) {
            complementarity_error += std::abs(multipliers.lower_bounds[i] * (iterate.x[i] - problem.variables_bounds[i].lb));
        }
        if (problem.variables_bounds[i].ub < INFINITY) {
            complementarity_error += std::abs(multipliers.upper_bounds[i] * (iterate.x[i] - problem.variables_bounds[i].ub));
        }
    }
    /* constraints */
    iterate.compute_constraints(problem);
    for (int j = 0; j < problem.number_constraints; j++) {
        double multiplier_j = multipliers.constraints[j];
        if (-INFINITY < problem.constraints_bounds[j].lb && 0. < multiplier_j) {
            complementarity_error += std::abs(multiplier_j * (iterate.constraints[j] - problem.constraints_bounds[j].lb));
        }
        if (problem.constraints_bounds[j].ub < INFINITY && multiplier_j < 0.) {
            complementarity_error += std::abs(multiplier_j * (iterate.constraints[j] - problem.constraints_bounds[j].ub));
        }
    }
    return complementarity_error;
}