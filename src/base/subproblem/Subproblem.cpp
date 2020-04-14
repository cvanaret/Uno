#include "Subproblem.hpp"

Subproblem::Subproblem() : number_subproblems_solved(0) {
}

Subproblem::~Subproblem() {
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
    current_iterate.compute_objective_gradient(problem);
    for (std::pair<int, double> term : current_iterate.objective_gradient) {
        int variable_index = term.first;
        double derivative = term.second;
        rhs[variable_index] += problem.objective_sign*derivative;
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