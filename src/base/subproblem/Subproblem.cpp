#include "Subproblem.hpp"
#include "LinearSolverFactory.hpp"

Subproblem::Subproblem(std::string residual_norm, std::vector<Range>& variables_bounds, bool scale_residuals):
residual_norm(residual_norm),
subproblem_variables_bounds(variables_bounds), // register the original bounds
number_subproblems_solved(0), subproblem_definition_changed(false),
scale_residuals(scale_residuals) {
}

Subproblem::~Subproblem() {
}

/* compute least-square multipliers */
//if (0 < problem.number_constraints) {
//    first_iterate.compute_constraints_jacobian(problem);
//    first_iterate.multipliers.constraints = Subproblem::compute_least_square_multipliers(problem, first_iterate, multipliers.constraints, 1e4);
//}

void Subproblem::project_point_in_bounds(std::vector<double>& x, const std::vector<Range>& variables_bounds) {
    for (unsigned int i = 0; i < x.size(); i++) {
        if (x[i] < variables_bounds[i].lb) {
            x[i] = variables_bounds[i].lb;
        }
        else if (variables_bounds[i].ub < x[i]) {
            x[i] = variables_bounds[i].ub;
        }
    }
    return;
}

double Subproblem::project_strictly_variable_in_bounds(double variable_value, const Range& variable_bounds) {
    double k1 = 1e-2;
    double k2 = 1e-2;

    double perturbation_lb = std::min(k1 * std::max(1., std::abs(variable_bounds.lb)), k2 * (variable_bounds.ub - variable_bounds.lb));
    double perturbation_ub = std::min(k1 * std::max(1., std::abs(variable_bounds.ub)), k2 * (variable_bounds.ub - variable_bounds.lb));
    variable_value = std::max(variable_value, variable_bounds.lb + perturbation_lb);
    variable_value = std::min(variable_value, variable_bounds.ub - perturbation_ub);
    return variable_value;
}

std::vector<Range> Subproblem::generate_constraints_bounds(const Problem& problem, const std::vector<double>& current_constraints) {
    std::vector<Range> constraints_bounds(problem.number_constraints);
    for (int j = 0; j < problem.number_constraints; j++) {
        double lb = problem.constraint_bounds[j].lb - current_constraints[j];
        double ub = problem.constraint_bounds[j].ub - current_constraints[j];
        constraints_bounds[j] = {lb, ub};
    }
    return constraints_bounds;
}

std::vector<double> Subproblem::compute_least_square_multipliers(const Problem& problem, Iterate& current_iterate, const std::vector<double>& default_multipliers, double multipliers_max_size) {
    std::unique_ptr<LinearSolver> linear_solver = LinearSolverFactory::create("MA57");
    return Subproblem::compute_least_square_multipliers(problem, current_iterate, default_multipliers, *linear_solver, multipliers_max_size);
}

std::vector<double> Subproblem::compute_least_square_multipliers(const Problem& problem, Iterate& current_iterate, const std::vector<double>& default_multipliers, LinearSolver& solver, double multipliers_max_size) {
    current_iterate.compute_objective_gradient(problem);
    current_iterate.compute_constraints_jacobian(problem);

    /******************************/
    /* build the symmetric matrix */
    /******************************/
    COOMatrix matrix(current_iterate.x.size() + problem.number_constraints, 1);

    /* identity blocks */
    for (unsigned int i = 0; i < current_iterate.x.size(); i++) {
        matrix.insert(1., i, i);
    }
    /* Jacobian of general constraints */
    for (int j = 0; j < problem.number_constraints; j++) {
        for (std::pair<const int, double>& term: current_iterate.constraints_jacobian[j]) {
            int variable_index = term.first;
            double derivative = term.second;
            matrix.insert(derivative, variable_index, current_iterate.x.size() + j);
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
    for (std::pair<int, double> term: current_iterate.objective_gradient) {
        int i = term.first;
        double derivative = term.second;
        rhs[i] += problem.objective_sign*derivative;
    }
    /* variable bound constraints */
    for (unsigned int i = 0; i < current_iterate.x.size(); i++) {
        rhs[i] -= current_iterate.multipliers.lower_bounds[i];
        rhs[i] -= current_iterate.multipliers.upper_bounds[i];
    }
    DEBUG << "Multipliers RHS:\n";
    print_vector(DEBUG, rhs);
    
    solver.do_symbolic_factorization(matrix);
    solver.solve(rhs);
    DEBUG << "Solution: ";
    std::vector<double>& solution = rhs;
    print_vector(DEBUG, solution);

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


double Subproblem::compute_KKT_error(const Problem& problem, Iterate& iterate, double objective_mutiplier) const {
    std::vector<double> lagrangian_gradient = iterate.lagrangian_gradient(problem, objective_mutiplier, iterate.multipliers);
    return norm(lagrangian_gradient, this->residual_norm);
}

/* complementary slackness error. Use abs/1e-8 to safeguard */
double Subproblem::compute_complementarity_error_(const Problem& problem, Iterate& iterate, const Multipliers& multipliers) const {
    double complementarity_error = 0.;
    /* bound constraints */
    for (unsigned int i = 0; i < iterate.x.size(); i++) {
        if (-INFINITY < this->subproblem_variables_bounds[i].lb) {
            complementarity_error += std::abs(multipliers.lower_bounds[i] * (iterate.x[i] - this->subproblem_variables_bounds[i].lb));
        }
        if (this->subproblem_variables_bounds[i].ub < INFINITY) {
            complementarity_error += std::abs(multipliers.upper_bounds[i] * (iterate.x[i] - this->subproblem_variables_bounds[i].ub));
        }
    }
    /* constraints */
    iterate.compute_constraints(problem);
    for (int j = 0; j < problem.number_constraints; j++) {
        double multiplier_j = multipliers.constraints[j];
        if (-INFINITY < problem.constraint_bounds[j].lb && 0. < multiplier_j) {
            complementarity_error += std::abs(multiplier_j * (iterate.constraints[j] - problem.constraint_bounds[j].lb));
        }
        if (problem.constraint_bounds[j].ub < INFINITY && multiplier_j < 0.) {
            complementarity_error += std::abs(multiplier_j * (iterate.constraints[j] - problem.constraint_bounds[j].ub));
        }
    }
    return complementarity_error;
}

void Subproblem::compute_residuals(const Problem& problem, Iterate& iterate, const Multipliers& multipliers, double objective_multiplier) const {
    iterate.compute_constraints(problem);
    iterate.residuals.constraints = problem.compute_constraint_residual(iterate.constraints, this->residual_norm);
    iterate.residuals.KKT = Subproblem::compute_KKT_error(problem, iterate, objective_multiplier);
    iterate.residuals.FJ = Subproblem::compute_KKT_error(problem, iterate, 0.);
    iterate.residuals.complementarity = this->compute_complementarity_error_(problem, iterate, multipliers);
    if (this->scale_residuals) {
        // TODO scale the residuals
    }
    return;
}
