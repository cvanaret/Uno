#include <cmath>
#include "InteriorPoint.hpp"
#include "Argonot.hpp"

InteriorPoint::InteriorPoint() : Subproblem(), mu(0.1), inertia_hessian(0.), inertia_hessian_last(0.), inertia_constraints(0.),
tau_min(0.99), default_multiplier(1.), k_sigma(1e10), smax(100.), k_mu(0.2), theta_mu(1.5), k_epsilon(10.), multipliers_max_size(1e3), iteration(0) {
}

Iterate InteriorPoint::initialize(Problem& problem, std::vector<double>& x, Multipliers& /* multipliers */, int /*number_variables*/, int /*number_constraints*/, bool use_trust_region) {
    int number_variables = problem.number_variables + problem.inequality_constraints.size();
    Multipliers multipliers(number_variables, problem.number_constraints);
    
    /* identify the variables' bounds */
    for (int i = 0; i < problem.number_variables; i++) {
        if (use_trust_region || (problem.variable_status[i] == BOUNDED_LOWER || problem.variable_status[i] == BOUNDED_BOTH_SIDES)) {
            this->lower_bounded_variables.push_back(i);
            multipliers.lower_bounds[i] = this->default_multiplier; // positive multiplier
        }
        if (use_trust_region || (problem.variable_status[i] == BOUNDED_UPPER || problem.variable_status[i] == BOUNDED_BOTH_SIDES)) {
            this->upper_bounded_variables.push_back(i);
            multipliers.upper_bounds[i] = -this->default_multiplier; // negative multiplier
        }
    }
    /* identify the inequality constraint slacks */
    for (std::pair<const int, int>& element: problem.inequality_constraints) {
        int j = element.first;
        int slack_index = problem.number_variables + element.second;
        if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
            this->lower_bounded_slacks[slack_index] = j;
            multipliers.lower_bounds[slack_index] = this->default_multiplier; // positive multiplier
        }
        if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
            this->upper_bounded_slacks[slack_index] = j;
            multipliers.upper_bounds[slack_index] = -this->default_multiplier; // negative multiplier
        }
    }
    
    /* make the initial point strictly feasible */
    std::vector<double> projected_x(number_variables);
    for (int i = 0; i < problem.number_variables; i++) {
        projected_x[i] = this->project_variable_in_bounds(x[i], problem.variables_bounds[i]);
    }

    /* generate the first iterate */
    Iterate first_iterate(problem, projected_x, multipliers);
    
    /* initialize the slacks */
    for (std::pair<const int, int>& element : problem.inequality_constraints) {
        int j = element.first;
        int slack_index = problem.number_variables + element.second;
        double slack_value = this->project_variable_in_bounds(first_iterate.constraints[j], problem.constraints_bounds[j]);
        first_iterate.x[slack_index] = slack_value;
    }
    /* compute least-square multipliers */
    first_iterate.multipliers.constraints = this->estimate_initial_multipliers(problem, first_iterate);
    
    DEBUG << problem.inequality_constraints.size() << " slacks\n";
    DEBUG << first_iterate.multipliers.lower_bounds.size() << " bound multipliers\n";
    DEBUG << first_iterate.multipliers.constraints.size() << " constraint multipliers\n";
    DEBUG << "variable lb: "; print_vector(DEBUG, this->lower_bounded_variables);
    DEBUG << "variable ub: "; print_vector(DEBUG, this->upper_bounded_variables);
    
    /* compute the optimality and feasibility measures of the initial point */
    this->compute_measures(problem, first_iterate);
    return first_iterate;
}

double InteriorPoint::compute_KKT_error_scaling(Iterate& current_iterate) {
    /* KKT error */
    double norm_1_constraint_multipliers = norm_1(current_iterate.multipliers.constraints);
    double norm_1_bound_multipliers = norm_1(current_iterate.multipliers.lower_bounds) + norm_1(current_iterate.multipliers.upper_bounds);
    double sd = std::max(this->smax, (norm_1_constraint_multipliers + norm_1_bound_multipliers)/(current_iterate.x.size() + current_iterate.multipliers.constraints.size())) / this->smax;
    return sd;
}

/* reduced primal-dual approach */
SubproblemSolution InteriorPoint::compute_optimality_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds) {
    DEBUG << "\nCurrent iterate: "; DEBUG << current_iterate;
    
    /* scaled error terms */
    double sd = this->compute_KKT_error_scaling(current_iterate);
    double KKTerror = Argonot::compute_KKT_error(problem, current_iterate, 1., "inf")/sd;
    double central_complementarity_error = this->compute_central_complementarity_error(problem, current_iterate, this->mu);
    DEBUG << "IPM error (KKT: " << KKTerror << ", cmpl: " << central_complementarity_error << ", feas: " << current_iterate.residual << ")\n";
    
    /* update of the barrier problem */
    double error = std::max(KKTerror, std::max(central_complementarity_error, current_iterate.residual));
    if (error <= this->k_epsilon*this->mu) {
        // TODO pass tolerance
        double tolerance = 1e-8;
        this->mu = std::max(tolerance/10., std::min(this->k_mu * this->mu, std::pow(this->mu, this->theta_mu)));
        DEBUG << "IPM: mu updated to " << this->mu << "\n";
    }
    DEBUG << "mu is " << this->mu << "\n";
    this->iteration++;

    /* compute first- and second-order information */
    current_iterate.compute_hessian(problem, problem.objective_sign, current_iterate.multipliers.constraints);
    current_iterate.compute_objective_gradient(problem);
    current_iterate.compute_constraints_jacobian(problem);

    /************************/
    /* solve the KKT system */
    /************************/
    /* KKT matrix */
    COOMatrix kkt_matrix = this->generate_kkt_matrix(problem, current_iterate, variables_bounds);

    /* right-hand side */
    std::vector<double> rhs = this->generate_kkt_rhs(problem, current_iterate, variables_bounds);

    /* compute the solution (Δx, -Δλ) */
    std::vector<double> solution_IPM = this->solver.solve(kkt_matrix, rhs, this->factorization_data);
    this->number_subproblems_solved++;

    /* generate IPM direction */
    SubproblemSolution solution = this->generate_direction(problem, current_iterate, solution_IPM, variables_bounds);
    solution.status = OPTIMAL;
    solution.norm = norm_inf(solution.x, problem.number_variables);

    /* evaluate the barrier objective */
    solution.objective = this->evaluate_local_model(problem, current_iterate, solution.x);
    return solution;
}

SubproblemSolution InteriorPoint::generate_direction(Problem& problem, Iterate& current_iterate, std::vector<double>& solution_IPM, std::vector<Range>& variables_bounds) {
    int number_variables = problem.number_variables + problem.inequality_constraints.size();
    
    /* retrieve +Δλ (Nocedal p590) */
    for (int j = 0; j < problem.number_constraints; j++) {
        int multiplier_index = number_variables + j;
        solution_IPM[multiplier_index] = -solution_IPM[multiplier_index];
    }

    /* compute bound multiplier displacements Δz */
    std::vector<double> lower_delta_z = this->compute_lower_bound_multiplier_displacements(problem, current_iterate, solution_IPM, variables_bounds);
    std::vector<double> upper_delta_z = this->compute_upper_bound_multiplier_displacements(problem, current_iterate, solution_IPM, variables_bounds);

    DEBUG << "MA57 solution:\n";
    DEBUG << "Δx: "; print_vector(DEBUG, solution_IPM, 0, problem.number_variables);
    DEBUG << "Δs: "; print_vector(DEBUG, solution_IPM, problem.number_variables, problem.inequality_constraints.size());
    DEBUG << "Δz_L: "; print_vector(DEBUG, lower_delta_z);
    DEBUG << "Δz_U: "; print_vector(DEBUG, upper_delta_z);
    DEBUG << "Δλ: "; print_vector(DEBUG, solution_IPM, number_variables, problem.number_constraints);

    /* "fraction to boundary" rule for variables and bound multipliers */
    std::vector<double> trial_x(current_iterate.x.size());
    Multipliers trial_multipliers(current_iterate.x.size(), current_iterate.constraints.size());
    double tau = std::max(this->tau_min, 1. - this->mu);
    // scale primal variables and constraints multipliers
    double primal_length = this->compute_primal_length(problem, current_iterate, solution_IPM, tau, variables_bounds);
    for (int i = 0; i < number_variables; i++) {
        trial_x[i] = primal_length * solution_IPM[i];
    }
    for (int j = 0; j < problem.number_constraints; j++) {
        trial_multipliers.constraints[j] = current_iterate.multipliers.constraints[j] + primal_length * solution_IPM[number_variables + j];
    }
    
    // scale dual variables
    double dual_length = this->compute_dual_length(current_iterate, tau, lower_delta_z, upper_delta_z);
    for (unsigned int i = 0; i < current_iterate.multipliers.lower_bounds.size(); i++) {
        trial_multipliers.lower_bounds[i] = current_iterate.multipliers.lower_bounds[i] + dual_length * lower_delta_z[i];
        trial_multipliers.upper_bounds[i] = current_iterate.multipliers.upper_bounds[i] + dual_length * upper_delta_z[i];
        // TODO rescale the multipliers (IPOPT paper p6)
        /*
        double scaling_ub = this->k_sigma * this->mu / trial_x[i];
        double scaling_lb = this->mu / (this->k_sigma * trial_x[i]);
        trial_lower_bound_multipliers[i] = std::max(std::min(trial_lower_bound_multipliers[i], scaling_ub), scaling_lb);
        trial_upper_bound_multipliers[i] = std::max(std::min(trial_upper_bound_multipliers[i], scaling_ub), scaling_lb);
         */
    }

    DEBUG << "primal length = " << primal_length << "\n";
    DEBUG << "dual length = " << dual_length << "\n\n";

    ActiveSet active_set;
    ConstraintPartition constraint_partition;
    SubproblemSolution solution(trial_x, trial_multipliers, active_set, constraint_partition);
    return solution;
}

std::vector<double> InteriorPoint::estimate_initial_multipliers(Problem& problem, Iterate& current_iterate) {
    /******************************/
    /* build the symmetric matrix */
    /******************************/
    COOMatrix matrix(problem.number_variables + problem.number_constraints, 0);
    
    /* identity blocks */
    for (int i = 0; i < problem.number_variables; i++) {
        matrix.add_term(1., i, i);
    }
    for (unsigned int slack_index = 0; slack_index < problem.inequality_constraints.size(); slack_index++) {
        matrix.add_term(-1., problem.number_variables + slack_index, problem.number_variables + slack_index);
    }

    /* Jacobian of general constraints */
    current_iterate.compute_constraints_jacobian(problem);
    for (int j = 0; j < problem.number_constraints; j++) {
        for (std::pair<const int, double>& term : current_iterate.constraints_jacobian[j]) {
            int variable_index = term.first;
            double derivative = term.second;
            matrix.add_term(derivative, variable_index, problem.number_variables + j);
        }
    }
    DEBUG << "Multipliers matrix:\n";
    for (unsigned int k = 0; k < matrix.matrix.size(); k++) {
        DEBUG << "m(" << matrix.row_indices[k] << ", " << matrix.column_indices[k] << ") = " << matrix.matrix[k] << "\n";
    }
    
    /********************************/
    /* generate the right-hand side */
    /********************************/
    std::vector<double> rhs(problem.number_variables + problem.number_constraints);

    /* objective gradient */
    current_iterate.compute_objective_gradient(problem);
    for (std::pair<int, double> term : current_iterate.objective_gradient) {
        int variable_index = term.first;
        double derivative = term.second;
        rhs[variable_index] += problem.objective_sign*derivative;
    }
    /* variable bound constraints */
    for (int i: this->lower_bounded_variables) {
        rhs[i] -= current_iterate.multipliers.lower_bounds[i];
    }
    for (int i: this->upper_bounded_variables) {
        rhs[i] -= current_iterate.multipliers.upper_bounds[i];
    }
    
    /* slack bound constraints */
    for (std::pair<const int, int>& element: this->lower_bounded_slacks) {
        int slack_index = element.first;
        int j = element.second;
        rhs[problem.number_variables + j] -= current_iterate.multipliers.lower_bounds[slack_index];
    }
    for (std::pair<const int, int>& element: this->upper_bounded_slacks) {
        int slack_index = element.first;
        int j = element.second;
        rhs[problem.number_variables + j] -= current_iterate.multipliers.upper_bounds[slack_index];
    }
    
    DEBUG << "Multipliers RHS:\n"; print_vector(DEBUG, rhs);
    
    MA57Data factorization_data = this->solver.factorize(matrix);
    std::vector<double> solution = this->solver.solve(matrix, rhs, factorization_data);
    
    /* retrieve multipliers */
    std::vector<double> multipliers(problem.number_constraints);
    for (int j = 0; j < problem.number_constraints; j++) {
        multipliers[j] = solution[problem.number_variables + j];
    }
    // if multipliers too big, discard them
    if (norm_inf(multipliers) > this->multipliers_max_size) {
        std::vector<double> empty_multipliers(problem.number_constraints);
        return empty_multipliers;
    }
    return multipliers;
}

double InteriorPoint::project_variable_in_bounds(double variable_value, Range& variable_bounds) {
    double k1 = 1e-2;
    double k2 = 1e-2;

    double perturbation_lb = std::min(k1 * std::max(1., std::abs(variable_bounds.lb)), k2 * (variable_bounds.ub - variable_bounds.lb));
    double perturbation_ub = std::min(k1 * std::max(1., std::abs(variable_bounds.ub)), k2 * (variable_bounds.ub - variable_bounds.lb));
    variable_value = std::max(variable_value, variable_bounds.lb + perturbation_lb);
    variable_value = std::min(variable_value, variable_bounds.ub - perturbation_ub);
    return variable_value;
}

double InteriorPoint::compute_primal_length(Problem& problem, Iterate& current_iterate, std::vector<double>& ipm_solution, double tau, std::vector<Range>& variables_bounds) {
    double primal_length = 1.;
    for (int i : this->lower_bounded_variables) {
        double trial_alpha_xi = -tau * (current_iterate.x[i] - variables_bounds[i].lb) / ipm_solution[i];
        if (0 < trial_alpha_xi && trial_alpha_xi <= 1.) {
            primal_length = std::min(primal_length, trial_alpha_xi);
        }
    }
    for (int i : this->upper_bounded_variables) {
        double trial_alpha_xi = -tau * (current_iterate.x[i] - variables_bounds[i].ub) / ipm_solution[i];
        if (0 < trial_alpha_xi && trial_alpha_xi <= 1.) {
            primal_length = std::min(primal_length, trial_alpha_xi);
        }
    }
    for (std::pair<const int, int>& element: this->lower_bounded_slacks) {
        int slack_index = element.first;
        int j = element.second;
        double trial_primal_length_sj = -tau * (current_iterate.x[slack_index] - problem.constraints_bounds[j].lb) / ipm_solution[slack_index];
        if (0 < trial_primal_length_sj && trial_primal_length_sj <= 1.) {
            primal_length = std::min(primal_length, trial_primal_length_sj);
        }
    }
    for (std::pair<const int, int>& element: this->upper_bounded_slacks) {
        int slack_index = element.first;
        int j = element.second;
        double trial_primal_length_sj = -tau * (current_iterate.x[slack_index] - problem.constraints_bounds[j].ub) / ipm_solution[slack_index];
        if (0 < trial_primal_length_sj && trial_primal_length_sj <= 1.) {
            primal_length = std::min(primal_length, trial_primal_length_sj);
        }
    }
    return primal_length;
}

double InteriorPoint::compute_dual_length(Iterate& current_iterate, double tau, std::vector<double>& lower_delta_z, std::vector<double>& upper_delta_z) {
    double dual_length = 1.;
    for (unsigned int current_multiplier = 0; current_multiplier < current_iterate.multipliers.lower_bounds.size(); current_multiplier++) {
        double trial_alpha_zj = -tau * current_iterate.multipliers.lower_bounds[current_multiplier] / lower_delta_z[current_multiplier];
        if (0 < trial_alpha_zj && trial_alpha_zj <= 1.) {
            dual_length = std::min(dual_length, trial_alpha_zj);
        }
        trial_alpha_zj = -tau * current_iterate.multipliers.upper_bounds[current_multiplier] / upper_delta_z[current_multiplier];
        if (0 < trial_alpha_zj && trial_alpha_zj <= 1.) {
            dual_length = std::min(dual_length, trial_alpha_zj);
        }
    }
    return dual_length;
}

COOMatrix InteriorPoint::generate_kkt_matrix(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds) {
    int number_variables = problem.number_variables + problem.inequality_constraints.size();
    
    /* compute the Lagrangian Hessian */
    current_iterate.compute_hessian(problem, problem.objective_sign, current_iterate.multipliers.constraints);
    COOMatrix kkt_matrix = current_iterate.hessian.to_COO();
    kkt_matrix.size = number_variables + problem.number_constraints;

    /* variable bound constraints */
    for (int i : this->lower_bounded_variables) {
        kkt_matrix.add_term(current_iterate.multipliers.lower_bounds[i] / (current_iterate.x[i] - variables_bounds[i].lb), i, i);
    }
    for (int i : this->upper_bounded_variables) {
        kkt_matrix.add_term(current_iterate.multipliers.upper_bounds[i] / (current_iterate.x[i] - variables_bounds[i].ub), i, i);
    }

    /* slack bound constraints */
    for (std::pair<const int, int>& element: this->lower_bounded_slacks) {
        int slack_index = element.first;
        int j = element.second;
        kkt_matrix.add_term(current_iterate.multipliers.lower_bounds[slack_index] / (current_iterate.x[slack_index] - problem.constraints_bounds[j].lb), slack_index, slack_index);
        kkt_matrix.add_term(-1., slack_index, number_variables + j);
    }
    for (std::pair<const int, int>& element: this->upper_bounded_slacks) {
        int slack_index = element.first;
        int j = element.second;
        kkt_matrix.add_term(current_iterate.multipliers.upper_bounds[slack_index] / (current_iterate.x[slack_index] - problem.constraints_bounds[j].ub), slack_index, slack_index);
        kkt_matrix.add_term(-1., slack_index, number_variables + j);
    }

    /* Jacobian of general constraints */
    for (int j = 0; j < problem.number_constraints; j++) {
        for (std::pair<int, double> term : current_iterate.constraints_jacobian[j]) {
            int variable_index = term.first;
            double derivative = term.second;
            kkt_matrix.add_term(derivative, variable_index, number_variables + j);
        }
    }

    /* inertia correction */
    this->modify_inertia(kkt_matrix, current_iterate.x.size(), problem.number_constraints);
    
    DEBUG << "KKT matrix:\n";
    for (unsigned int k = 0; k < kkt_matrix.matrix.size(); k++) {
        DEBUG << "m(" << kkt_matrix.row_indices[k] << ", " << kkt_matrix.column_indices[k] << ") = " << kkt_matrix.matrix[k] << "\n";
    }
    
    return kkt_matrix;
}

void InteriorPoint::modify_inertia(COOMatrix& kkt_matrix, int number_variables, int number_constraints) {
    this->inertia_hessian = 0.;
    this->inertia_constraints = 0.;
    DEBUG << "Testing factorization with inertia term " << this->inertia_hessian << "\n";
    this->factorization_data = this->solver.factorize(kkt_matrix);

    bool good_inertia = false;
    if (!this->solver.matrix_is_singular() && this->solver.number_negative_eigenvalues() == number_constraints) {
        DEBUG << "Factorization was a success\n";
        good_inertia = true;
    }
    else {
        // inertia term for constraints
        if (this->solver.matrix_is_singular()) {
            DEBUG << "Matrix is singular\n";
            this->inertia_constraints = 1e-8 * std::pow(this->mu, 0.25);
        }
        else {
            this->inertia_constraints = 0.;
        }
        // inertia term for Hessian
        if (this->inertia_hessian_last == 0.) {
            this->inertia_hessian = 1e-4;
        } else {
            this->inertia_hessian = std::max(1e-20, this->inertia_hessian_last / 3.);
        }
    }

    int current_matrix_size = kkt_matrix.matrix.size();
    if (!good_inertia) {
        for (int i = 0; i < number_variables; i++) {
            kkt_matrix.add_term(this->inertia_hessian, i, i);
        }
        for (int j = number_variables; j < number_variables + number_constraints; j++) {
            kkt_matrix.add_term(-this->inertia_constraints, j, j);
        }
    }

    while (!good_inertia) {
        DEBUG << "Testing factorization with inertia term " << this->inertia_hessian << "\n";
        this->factorization_data = this->solver.factorize(kkt_matrix);

        if (!this->solver.matrix_is_singular() && this->solver.number_negative_eigenvalues() == number_constraints) {
            good_inertia = true;
            DEBUG << "Factorization was a success\n";
            this->inertia_hessian_last = this->inertia_hessian;
        }
        else {
            if (this->inertia_hessian_last == 0.) {
                this->inertia_hessian *= 100.;
            }
            else {
                this->inertia_hessian *= 8.;
            }
            if (1e40 < this->inertia_hessian) {
                throw std::out_of_range("Switch to restoration phase");
            }
            else {
                for (int i = 0; i < number_variables; i++) {
                    kkt_matrix.matrix[current_matrix_size + i] = this->inertia_hessian;
                }
                for (int j = number_variables; j < number_variables + number_constraints; j++) {
                    kkt_matrix.matrix[current_matrix_size + j] = -this->inertia_constraints;
                }
            }
        }
    }
    return;
}

std::vector<double> InteriorPoint::generate_kkt_rhs(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds) {
    int number_variables = problem.number_variables + problem.inequality_constraints.size();
    
    /* generate the right-hand side */
    std::vector<double> rhs(number_variables + problem.number_constraints);
    
    /* objective gradient */
    for (std::pair<int, double> term : current_iterate.objective_gradient) {
        int variable_index = term.first;
        double derivative = term.second;
        rhs[variable_index] = -problem.objective_sign*derivative;
    }
    /* constraint gradients */
    for (int j = 0; j < problem.number_constraints; j++) {
        if (current_iterate.multipliers.constraints[j] != 0.) {
            for (std::pair<int, double> term : current_iterate.constraints_jacobian[j]) {
                int variable_index = term.first;
                double derivative = term.second;
                rhs[variable_index] += current_iterate.multipliers.constraints[j]*derivative;
            }
        }
    }

    /* variable bound constraints */
    for (int i : this->lower_bounded_variables) {
        rhs[i] += this->mu / (current_iterate.x[i] - variables_bounds[i].lb);
    }
    for (int i : this->upper_bounded_variables) {
        rhs[i] += this->mu / (current_iterate.x[i] - variables_bounds[i].ub);
    }

    /* slack bound constraints */
    for (std::pair<const int, int>& element: this->lower_bounded_slacks) {
        int slack_index = element.first;
        int j = element.second;
        rhs[slack_index] = this->mu / (current_iterate.x[slack_index] - problem.constraints_bounds[j].lb) - current_iterate.multipliers.constraints[j];
    }
    for (std::pair<const int, int>& element: this->upper_bounded_slacks) {
        int slack_index = element.first;
        int j = element.second;
        rhs[slack_index] = this->mu / (current_iterate.x[slack_index] - problem.constraints_bounds[j].ub) - current_iterate.multipliers.constraints[j];
    }
    
    /* constraint evaluations */
    int slack_index = 0;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] == EQUAL_BOUNDS) {
            // add the bound
            rhs[number_variables + j] = -(current_iterate.constraints[j] - problem.constraints_bounds[j].lb);
        }
        else {
            // add the slack
            rhs[number_variables + j] = -(current_iterate.constraints[j] - current_iterate.x[problem.number_variables + slack_index]);
            slack_index++;
        }
    }
    DEBUG << "RHS: "; print_vector(DEBUG, rhs);

    return rhs;
}

std::vector<double> InteriorPoint::compute_lower_bound_multiplier_displacements(Problem& problem, Iterate& current_iterate, std::vector<double>& solution, std::vector<Range>& variables_bounds) {
    std::vector<double> delta_z(current_iterate.multipliers.lower_bounds.size());
    for (int i : this->lower_bounded_variables) {
        delta_z[i] = this->mu / (current_iterate.x[i] - variables_bounds[i].lb) - current_iterate.multipliers.lower_bounds[i] - current_iterate.multipliers.lower_bounds[i] / (current_iterate.x[i] - variables_bounds[i].lb) * solution[i];
    }
    for (std::pair<const int, int>& element: this->lower_bounded_slacks) {
        int slack_index = element.first;
        int j = element.second;
        delta_z[slack_index] = this->mu / (current_iterate.x[slack_index] - problem.constraints_bounds[j].lb) - current_iterate.multipliers.lower_bounds[slack_index] - current_iterate.multipliers.lower_bounds[slack_index] / (current_iterate.x[slack_index] - problem.constraints_bounds[j].lb) * solution[slack_index];
    }
    return delta_z;
}

std::vector<double> InteriorPoint::compute_upper_bound_multiplier_displacements(Problem& problem, Iterate& current_iterate, std::vector<double>& solution, std::vector<Range>& variables_bounds) {
    std::vector<double> delta_z(current_iterate.multipliers.upper_bounds.size());
    for (int i : this->upper_bounded_variables) {
        delta_z[i] = this->mu / (current_iterate.x[i] - variables_bounds[i].ub) - current_iterate.multipliers.upper_bounds[i] - current_iterate.multipliers.upper_bounds[i] / (current_iterate.x[i] - variables_bounds[i].ub) * solution[i];
    }
    for (std::pair<const int, int>& element: this->upper_bounded_slacks) {
        int slack_index = element.first;
        int j = element.second;
        delta_z[slack_index] = this->mu / (current_iterate.x[slack_index] - problem.constraints_bounds[j].ub) - current_iterate.multipliers.upper_bounds[slack_index] - current_iterate.multipliers.upper_bounds[slack_index] / (current_iterate.x[slack_index] - problem.constraints_bounds[j].ub) * solution[slack_index];
    }
    return delta_z;
}

void InteriorPoint::compute_measures(Problem& problem, Iterate& iterate) {
    /* evaluate constraints with slacks */
    iterate.feasibility_measure = this->constraint_violation(problem, iterate);
    /* compute barrier objective */
    iterate.optimality_measure = this->barrier_objective(problem, iterate);
    return;
}

double InteriorPoint::constraint_violation(Problem& problem, Iterate& iterate) {
    double constraint_violation = 0.;
    int slack_index = problem.number_variables;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] == EQUAL_BOUNDS) {
            double constraint_value = iterate.constraints[j] - problem.constraints_bounds[j].lb;
            constraint_violation += std::abs(constraint_value);
        }
        else {
            double constraint_value = iterate.constraints[j] - iterate.x[slack_index];
            constraint_violation += std::abs(constraint_value);
            slack_index++;
        }
    }
    return constraint_violation;
}

double InteriorPoint::barrier_objective(Problem& problem, Iterate& iterate) {
    /* original objective */
    double objective = iterate.objective;
    
    /* variable bound constraints */
    for (int i : this->lower_bounded_variables) {
        objective -= this->mu * std::log(iterate.x[i] - problem.variables_bounds[i].lb);
    }
    for (int i : this->upper_bounded_variables) {
        objective -= this->mu * std::log(problem.variables_bounds[i].ub - iterate.x[i]);
    }
    
    /* slack bound constraints */
    for (std::pair<const int, int>& element: this->lower_bounded_slacks) {
        int slack_index = element.first;
        int j = element.second;
        objective -= this->mu * std::log(iterate.x[slack_index] - problem.constraints_bounds[j].lb);
    }
    for (std::pair<const int, int>& element: this->upper_bounded_slacks) {
        int slack_index = element.first;
        int j = element.second;
        objective -= this->mu * std::log(problem.constraints_bounds[j].ub - iterate.x[slack_index]);
    }
    return objective;
}

std::map<int, double> InteriorPoint::barrier_function_gradient(Problem& problem, Iterate& current_iterate) {
    current_iterate.compute_objective_gradient(problem);
    /* objective gradient */
    std::map<int, double> gradient(current_iterate.objective_gradient);
    
    /* barrier terms for variables */
    for (int i : this->lower_bounded_variables) {
        gradient[i] -= this->mu / (current_iterate.x[i] - problem.variables_bounds[i].lb);
    }
    for (int i : this->upper_bounded_variables) {
        gradient[i] -= this->mu / (current_iterate.x[i] - problem.variables_bounds[i].ub);
    }
    
    /* barrier terms for slacks */
    for (std::pair<const int, int>& element: this->lower_bounded_slacks) {
        int slack_index = element.first;
        int j = element.second;
        gradient[slack_index] -= this->mu / (current_iterate.x[slack_index] - problem.constraints_bounds[j].lb);
    }
    for (std::pair<const int, int>& element: this->upper_bounded_slacks) {
        int slack_index = element.first;
        int j = element.second;
        gradient[slack_index] -= this->mu / (current_iterate.x[slack_index] - problem.constraints_bounds[j].ub);
    }
    return gradient;
}

double InteriorPoint::evaluate_local_model(Problem& problem, Iterate& current_iterate, std::vector<double>& solution) {
    std::map<int, double> gradient = this->barrier_function_gradient(problem, current_iterate);
    double subproblem_objective = dot(solution, gradient);
    return subproblem_objective;
}

double InteriorPoint::compute_predicted_reduction(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution, double step_length) {
    double subproblem_objective = this->evaluate_local_model(problem, current_iterate, solution.x);
    return -step_length*subproblem_objective;
}

SubproblemSolution InteriorPoint::compute_infeasibility_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds, SubproblemSolution& phase_II_solution) {
    throw std::out_of_range("ENTERING IPM.compute_infeasibility_step, no implementation provided!");
}

SubproblemSolution InteriorPoint::compute_l1_penalty_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds, double penalty_parameter, PenaltyDimensions penalty_dimensions) {
    throw std::out_of_range("ENTERING IPM.compute_l1_penalty_step, no implementation provided!");
}

bool InteriorPoint::phase_1_required(SubproblemSolution& solution) {
    return solution.status == INFEASIBLE;
}

double InteriorPoint::compute_central_complementarity_error(const Problem& problem, Iterate& iterate, double mu) {
    double complementarity_error = 0.;
    
    /* variable bound constraints */
    for (int i = 0; i < problem.number_variables; i++) {
        if (-INFINITY < problem.variables_bounds[i].lb) {
            complementarity_error = std::max(complementarity_error, std::abs(iterate.multipliers.lower_bounds[i] * (iterate.x[i] - problem.variables_bounds[i].lb) - mu));
        }
        if (problem.variables_bounds[i].ub < INFINITY) {
            complementarity_error = std::max(complementarity_error, std::abs(iterate.multipliers.upper_bounds[i] * (iterate.x[i] - problem.variables_bounds[i].ub) - mu));
        }
    }
    
    /* slack bound constraints */
    for (std::pair<const int, int>& element: this->lower_bounded_slacks) {
        int slack_index = element.first;
        int j = element.second;
        if (-INFINITY < problem.constraints_bounds[j].lb) {
            complementarity_error = std::max(complementarity_error, std::abs(iterate.multipliers.lower_bounds[slack_index] * (iterate.x[slack_index] - problem.constraints_bounds[j].lb) - mu));
        }
    }
    for (std::pair<const int, int>& element: this->upper_bounded_slacks) {
        int slack_index = element.first;
        int j = element.second;
        if (problem.constraints_bounds[j].ub < INFINITY) {
            complementarity_error = std::max(complementarity_error, std::abs(iterate.multipliers.upper_bounds[slack_index] * (iterate.x[slack_index] - problem.constraints_bounds[j].ub) - mu));
        }
    }
    
    /* scaling */
    double sc = std::max(this->smax, (norm_1(iterate.multipliers.lower_bounds) + norm_1(iterate.multipliers.upper_bounds))/iterate.x.size()) / this->smax;
    complementarity_error /= sc;
    return complementarity_error;
}