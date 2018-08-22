#include <cmath>
#include "InteriorPoint.hpp"

InteriorPoint::InteriorPoint() : Subproblem(), mu(0.1), inertia_hessian(0.), inertia_hessian_last(0.), inertia_constraints(0.),
tau_min(0.99), default_multiplier(1.), k_sigma(1e10), iteration(0) {
}

void InteriorPoint::initialize(Problem& problem, Iterate& first_iterate, int /*number_variables*/, int /*number_constraints*/, bool use_trust_region) {
    /* identify the variables' bounds */
    for (int i = 0; i < problem.number_variables; i++) {
        if (use_trust_region || (problem.variable_status[i] == BOUNDED_LOWER || problem.variable_status[i] == BOUNDED_BOTH_SIDES)) {
            this->lower_bounded_variables.push_back(i);
        }
        if (use_trust_region || (problem.variable_status[i] == BOUNDED_UPPER || problem.variable_status[i] == BOUNDED_BOTH_SIDES)) {
            this->upper_bounded_variables.push_back(i);
        }
    }
    /* identify the inequality constraint slacks */
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] != EQUAL_BOUNDS) {
            this->slacked_constraints.push_back(j);
            if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
                this->lower_bounded_slacks.push_back(j);
            }
            if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
                this->upper_bounded_slacks.push_back(j);
            }
        }
    }

    /* initialize the slacks and multipliers */
    std::vector<double> bound_multipliers;
    for (int i : this->lower_bounded_variables) {
        bound_multipliers.push_back(this->default_multiplier); // positive multiplier
    }
    for (int i : this->upper_bounded_variables) {
        bound_multipliers.push_back(-this->default_multiplier); // negative multiplier
    }
    for (int j : this->lower_bounded_slacks) {
        bound_multipliers.push_back(this->default_multiplier); // positive multiplier
        double slack_value = this->project_variable_in_bounds(first_iterate.constraints[j], problem.constraint_lb[j], problem.constraint_ub[j]);
        first_iterate.x.push_back(slack_value);
    }
    for (int j : this->upper_bounded_slacks) {
        bound_multipliers.push_back(-this->default_multiplier); // negative multiplier
        double slack_value = this->project_variable_in_bounds(first_iterate.constraints[j], problem.constraint_lb[j], problem.constraint_ub[j]);
        first_iterate.x.push_back(slack_value);
    }
    first_iterate.bound_multipliers = bound_multipliers;

    /* compute first-order information */
    first_iterate.compute_objective_gradient(problem);
    first_iterate.compute_constraints_jacobian(problem);
    /* compute least-square multipliers */
    first_iterate.constraint_multipliers = this->estimate_initial_multipliers(problem, first_iterate);

    std::cout << this->slacked_constraints.size() << " slacks\n";
    std::cout << first_iterate.bound_multipliers.size() << " bound multipliers\n";
    std::cout << first_iterate.constraint_multipliers.size() << " constraint multipliers\n";
    std::cout << "variable lb: ";
    print_vector(std::cout, this->lower_bounded_variables);
    std::cout << "variable ub: ";
    print_vector(std::cout, this->upper_bounded_variables);
    std::cout << "slacks: ";
    print_vector(std::cout, this->slacked_constraints);
    std::cout << "slack lb: ";
    print_vector(std::cout, this->lower_bounded_slacks);
    std::cout << "slack ub: ";
    print_vector(std::cout, this->upper_bounded_slacks);

    /* compute the optimality and feasibility measures of the initial point */
    this->compute_measures(problem, first_iterate);

    return;
}

/* reduced primal-dual approach */
LocalSolution InteriorPoint::compute_optimality_step(Problem& problem, Iterate& current_iterate, double radius) {
    /* build trust region */
    std::vector<double> variable_lb(problem.number_variables);
    std::vector<double> variable_ub(problem.number_variables);
    for (int i = 0; i < problem.number_variables; i++) {
        variable_lb[i] = std::max(current_iterate.x[i] - radius, problem.variable_lb[i]);
        variable_ub[i] = std::min(current_iterate.x[i] + radius, problem.variable_ub[i]);
    }

    /* update barrier parameter */
    // double k_mu = 0.2;
    // double theta_mu = 1.5;
    // this->mu = std::max(1e-7, std::min(k_mu * this->mu, std::pow(this->mu, theta_mu)));
    // this->mu = this->update_barrier_parameter(problem, current_iterate);
    // std::cout << "NEW mu VALUE: " << this->mu << "\n";

    if (iteration == 1) {
        this->mu = 2.828427e-3;
    }
    else if (iteration == 2) {
        this->mu = 1.504241e-04;
    }
    else if (iteration == 4) {
        this->mu = 1.844914e-06;
    }
    else if (iteration == 6) {
        this->mu = 2.505904e-09;
    }

    this->iteration++;

    /* summary */
    //    for (int i = 0; i < problem.number_variables; i++) {
    //        std::cout << "x" << i << " in [" << variable_lb[i] << ", " << variable_ub[i] << "]\n";
    //    }
    //    std::cout << "\n";
    //    for (int j = 0; j < problem.number_constraints; j++) {
    //        std::cout << "c" << j << " in [" << problem.constraint_lb[j] << ", " << problem.constraint_ub[j] << "]\n";
    //    }
    //    std::cout << "\n";
    //    std::cout << "x/s is: ";
    //    print_vector(std::cout, current_iterate.x);
    //    std::cout << "λ is: ";
    //    print_vector(std::cout, current_iterate.constraint_multipliers);
    //    std::cout << "z is: ";
    //    print_vector(std::cout, current_iterate.bound_multipliers);
    //    std::cout << "Constraints: ";
    //    print_vector(std::cout, current_iterate.constraints);
    //    std::cout << "\n";

    std::cout << "mu is " << this->mu << "\n";

    /* compute first-order information */
    current_iterate.compute_objective_gradient(problem);
    current_iterate.compute_constraints_jacobian(problem);

    /*******************************/
    /* sparse symmetric KKT matrix */
    /*******************************/
    // create the KKT matrix
    COOMatrix kkt_matrix = this->generate_kkt_matrix(problem, current_iterate, variable_lb, variable_ub);

    /*******************/
    /* right-hand side */
    /*******************/
    std::vector<double> rhs = this->generate_rhs(problem, current_iterate, variable_lb, variable_ub);

    /************/
    /* solution */
    /************/
    /* compute the solution (Δx, -Δλ) */
    std::vector<double> ipm_solution = this->solver.solve(kkt_matrix, rhs, this->data);
    this->number_subproblems_solved++;

    /* retrieve +Δλ (Nocedal p590) */
    for (int j = 0; j < problem.number_constraints; j++) {
        int multiplier_index = problem.number_variables + this->slacked_constraints.size() + j;
        ipm_solution[multiplier_index] = -ipm_solution[multiplier_index];
    }

    /* compute bound multiplier displacements Δz */
    std::vector<double> delta_z = this->compute_bound_multiplier_displacements(problem, current_iterate, ipm_solution, variable_lb, variable_ub);

    std::cout << "MA57 solution:\n";
    std::cout << "Δx: ";
    print_vector(std::cout, ipm_solution, 0, problem.number_variables);
    std::cout << "Δs: ";
    print_vector(std::cout, ipm_solution, problem.number_variables, this->slacked_constraints.size());
    std::cout << "Δz: ";
    print_vector(std::cout, delta_z);
    std::cout << "Δλ: ";
    print_vector(std::cout, ipm_solution, problem.number_variables + this->slacked_constraints.size(), problem.number_constraints);

    /* "fraction to boundary" rule for variables and bound multipliers */
    double tau = std::max(this->tau_min, 1. - this->mu);
    double primal_length = this->compute_primal_length(problem, current_iterate, ipm_solution, tau, variable_lb, variable_ub);
    std::cout << "primal length = " << primal_length << "\n";
    double dual_length = this->compute_dual_length(current_iterate, tau, delta_z);
    std::cout << "dual length = " << dual_length << "\n\n";

    /* generate IPM direction */
    std::vector<double> trial_x(current_iterate.x.size());
    std::vector<double> trial_bound_multipliers(current_iterate.bound_multipliers.size());
    std::vector<double> trial_constraint_multipliers(current_iterate.constraint_multipliers.size());

    for (int i = 0; i < problem.number_variables; i++) {
        trial_x[i] = primal_length * ipm_solution[i];
    }
    for (unsigned int current_slack = 0; current_slack < this->slacked_constraints.size(); current_slack++) {
        trial_x[problem.number_variables + current_slack] = primal_length * ipm_solution[problem.number_variables + current_slack];
    }
    std::cout << "Scaled Δx/Δs: ";
    print_vector(std::cout, trial_x);

    for (unsigned int i = 0; i < current_iterate.bound_multipliers.size(); i++) {
        trial_bound_multipliers[i] = current_iterate.bound_multipliers[i] + dual_length * delta_z[i];
    }
    /* reset z */
    //    current_multiplier = 0;
    //    for (int i = 0; i < problem.number_variables; i++) {
    //        if (variable_status[i] == BOUNDED_LOWER || variable_status[i] == BOUNDED_BOTH_SIDES) {
    //            std::max(std::min(current_iterate.bound_multipliers[current_multiplier], this->k_sigma * this->mu / trial_primal[i]), this->mu / (this->k_sigma * trial_primal[i]));
    //            current_multiplier++;
    //        }
    //        if (variable_status[i] == BOUNDED_UPPER || variable_status[i] == BOUNDED_BOTH_SIDES) {
    //            std::max(std::min(current_iterate.bound_multipliers[current_multiplier], this->k_sigma * this->mu / trial_primal[i]), this->mu / (this->k_sigma * trial_primal[i]));
    //            current_multiplier++;
    //        }
    //    }
    // TODO reset z for s
    std::cout << "Scaled z: ";
    print_vector(std::cout, trial_bound_multipliers);

    for (int j = 0; j < problem.number_constraints; j++) {
        trial_constraint_multipliers[j] = primal_length * ipm_solution[problem.number_variables + this->slacked_constraints.size() + j];
    }
    std::cout << "Scaled Δλ: ";
    print_vector(std::cout, trial_constraint_multipliers);

    LocalSolution solution(trial_x, trial_bound_multipliers, trial_constraint_multipliers);
    solution.status = OPTIMAL;
    solution.norm = norm_inf(solution.x, problem.number_variables);
    /* compute */
    solution.objective = this->evaluate_local_model(problem, current_iterate, trial_x);
    solution.objective_terms.linear = solution.objective;
    solution.objective_terms.quadratic = 0.;
    
    return solution;
}

std::vector<double> InteriorPoint::estimate_initial_multipliers(Problem& problem, Iterate& current_iterate) {
    /* build the symmetric matrix */
    COOMatrix matrix(problem.number_variables + problem.number_constraints, 0);

    /* identity blocks */
    for (int i = 0; i < problem.number_variables; i++) {
        matrix.add_term(-1., i, i);
    }
    for (unsigned int current_slack = 0; current_slack < this->slacked_constraints.size(); current_slack++) {
        matrix.add_term(1., problem.number_variables + current_slack, problem.number_variables + current_slack);
    }
    /* Jacobian of general constraints */
    int current_column = problem.number_variables;
    for (int j = 0; j < problem.number_constraints; j++) {
        for (std::pair<int, double> term : current_iterate.constraints_jacobian[j]) {
            int variable_index = term.first;
            double derivative = term.second;
            matrix.add_term(derivative, variable_index, current_column);
        }
        current_column++;
    }

    DEBUG << "Matrix:\n";
    for (unsigned int k = 0; k < matrix.matrix.size(); k++) {
        DEBUG << "m(" << matrix.row_indices[k] << ", " << matrix.column_indices[k] << ") = " << matrix.matrix[k] << "\n";
    }

    /* generate the right-hand side */
    std::vector<double> rhs(problem.number_variables + problem.number_constraints);

    /* objective gradient */
    for (std::pair<int, double> term : current_iterate.objective_gradient) {
        int variable_index = term.first;
        double derivative = term.second;
        rhs[variable_index] += problem.objective_sign*derivative;
    }
    /* variable bound constraints */
    int current_multiplier = 0;
    for (int i : this->lower_bounded_variables) {
        rhs[i] -= current_iterate.bound_multipliers[current_multiplier];
        current_multiplier++;
    }
    for (int i : this->upper_bounded_variables) {
        rhs[i] -= current_iterate.bound_multipliers[current_multiplier];
        current_multiplier++;
    }

    /* slack bound constraints */
    for (int j : this->lower_bounded_slacks) {
        DEBUG << "rhs[" << (problem.number_variables + j) << "] += " << current_iterate.bound_multipliers[current_multiplier] << "\n";
        rhs[problem.number_variables + j] += current_iterate.bound_multipliers[current_multiplier];
        current_multiplier++;
    }
    for (int j : this->upper_bounded_slacks) {
        DEBUG << "rhs[" << (problem.number_variables + j) << "] += " << current_iterate.bound_multipliers[current_multiplier] << "\n";
        rhs[problem.number_variables + j] += current_iterate.bound_multipliers[current_multiplier];
        current_multiplier++;
    }

    MA57Data data = this->solver.factorize(matrix);
    std::vector<double> solution = this->solver.solve(matrix, rhs, data);

    /* retrieve multipliers */
    std::vector<double> multipliers(problem.number_constraints);
    for (int j = 0; j < problem.number_constraints; j++) {
        multipliers[j] = solution[problem.number_variables + j];
    }
    return multipliers;
}

double InteriorPoint::project_variable_in_bounds(double current_value, double lb, double ub) {
    double k1 = 1e-2;
    double k2 = 1e-2;

    double perturbation_lb = std::min(k1 * std::max(1., std::abs(lb)), k2 * (ub - lb));
    double perturbation_ub = std::min(k1 * std::max(1., std::abs(ub)), k2 * (ub - lb));
    current_value = std::max(current_value, lb + perturbation_lb);
    current_value = std::min(current_value, ub - perturbation_ub);
    return current_value;
}

double InteriorPoint::compute_primal_length(Problem& problem, Iterate& current_iterate, std::vector<double>& ipm_solution, double tau, std::vector<double> variable_lb, std::vector<double> variable_ub) {
    double primal_length = 1.;
    for (int i : this->lower_bounded_variables) {
        double trial_alpha_xi = -tau * (current_iterate.x[i] - variable_lb[i]) / ipm_solution[i];
        if (0 < trial_alpha_xi && trial_alpha_xi <= 1.) {
            primal_length = std::min(primal_length, trial_alpha_xi);
        }
    }
    for (int i : this->upper_bounded_variables) {
        double trial_alpha_xi = -tau * (current_iterate.x[i] - variable_ub[i]) / ipm_solution[i];
        if (0 < trial_alpha_xi && trial_alpha_xi <= 1.) {
            primal_length = std::min(primal_length, trial_alpha_xi);
        }
    }
    int current_slack = 0;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
            double trial_primal_lengthj = -tau * (current_iterate.x[problem.number_variables + current_slack] - problem.constraint_lb[j]) / ipm_solution[problem.number_variables + current_slack];
            if (0 < trial_primal_lengthj && trial_primal_lengthj <= 1.) {
                primal_length = std::min(primal_length, trial_primal_lengthj);
            }
        }
        if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
            double trial_primal_lengthj = -tau * (current_iterate.x[problem.number_variables + current_slack] - problem.constraint_ub[j]) / ipm_solution[problem.number_variables + current_slack];
            if (0 < trial_primal_lengthj && trial_primal_lengthj <= 1.) {
                primal_length = std::min(primal_length, trial_primal_lengthj);
            }
        }
        if (problem.constraint_status[j] != EQUAL_BOUNDS) {
            current_slack++;
        }
    }
    return primal_length;
}

double InteriorPoint::compute_dual_length(Iterate& current_iterate, double tau, std::vector<double>& delta_z) {
    double dual_length = 1.;
    for (unsigned int current_multiplier = 0; current_multiplier < current_iterate.bound_multipliers.size(); current_multiplier++) {
        double trial_alpha_zj = -tau * current_iterate.bound_multipliers[current_multiplier] / delta_z[current_multiplier];
        if (0 < trial_alpha_zj && trial_alpha_zj <= 1.) {
            dual_length = std::min(dual_length, trial_alpha_zj);
        }
    }
    return dual_length;
}

COOMatrix InteriorPoint::generate_kkt_matrix(Problem& problem, Iterate& current_iterate, std::vector<double>& variable_lb, std::vector<double>& variable_ub) {
    /* compute the Lagrangian Hessian */
    current_iterate.compute_hessian(problem, problem.objective_sign, current_iterate.constraint_multipliers);
    COOMatrix kkt_matrix = current_iterate.hessian.to_COO();
    kkt_matrix.size = problem.number_variables + this->slacked_constraints.size() + problem.number_constraints;

    /* variable bound constraints */
    int current_multiplier = 0;
    for (int i : this->lower_bounded_variables) {
        kkt_matrix.add_term(current_iterate.bound_multipliers[current_multiplier] / (current_iterate.x[i] - variable_lb[i]), i, i);
        current_multiplier++;
    }
    for (int i : this->upper_bounded_variables) {
        kkt_matrix.add_term(current_iterate.bound_multipliers[current_multiplier] / (current_iterate.x[i] - variable_ub[i]), i, i);
        current_multiplier++;
    }

    /* slack bound constraints */
    int current_column = problem.number_variables;
    int current_slack = 0;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
            kkt_matrix.add_term(current_iterate.bound_multipliers[current_multiplier] / (current_iterate.x[problem.number_variables + current_slack] - problem.constraint_lb[j]), current_column, current_column);
            kkt_matrix.add_term(-1., current_column, problem.number_variables + this->slacked_constraints.size() + j);
            current_multiplier++;
        }
        if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
            kkt_matrix.add_term(current_iterate.bound_multipliers[current_multiplier] / (current_iterate.x[problem.number_variables + current_slack] - problem.constraint_ub[j]), current_column, current_column);
            kkt_matrix.add_term(-1., current_column, problem.number_variables + this->slacked_constraints.size() + j);
            current_multiplier++;
        }
        if (problem.constraint_status[j] != EQUAL_BOUNDS) {
            current_slack++;
            current_column++;
        }
    }

    /* Jacobian of general constraints */
    for (int j = 0; j < problem.number_constraints; j++) {
        for (std::pair<int, double> term : current_iterate.constraints_jacobian[j]) {
            int variable_index = term.first;
            double derivative = term.second;
            kkt_matrix.add_term(derivative, variable_index, current_column);
        }
        current_column++;
    }

    /* Inertia correction */
    this->inertia_correction(problem, kkt_matrix);

    return kkt_matrix;
}

void InteriorPoint::inertia_correction(Problem& problem, COOMatrix& kkt_matrix) {
    this->inertia_hessian = 0.;
    this->inertia_constraints = 0.;
    std::cout << "Testing factorization with inertia term " << this->inertia_hessian << "\n";
    DEBUG << "KKT matrix:\n";
    for (unsigned int k = 0; k < kkt_matrix.matrix.size(); k++) {
        DEBUG << "m(" << kkt_matrix.row_indices[k] << ", " << kkt_matrix.column_indices[k] << ") = " << kkt_matrix.matrix[k] << "\n";
    }
    this->data = this->solver.factorize(kkt_matrix);

    bool good_inertia = false;
    if (!this->solver.matrix_is_singular() && this->solver.number_negative_eigenvalues() == problem.number_constraints) {
        std::cout << "Factorization was a success\n";
        good_inertia = true;
    }
    else {
        // inertia term for constraints
        if (this->solver.matrix_is_singular()) {
            std::cout << "Matrix is singular\n";
            this->inertia_constraints = 1e-8 * std::pow(this->mu, 0.25);
        }
        else {
            this->inertia_constraints = 0.;
        }
        // inertia term for Hessian
        if (this->inertia_hessian_last == 0.) {
            this->inertia_hessian = 1e-4;
        }
        else {
            this->inertia_hessian = std::max(1e-20, this->inertia_hessian_last / 3.);
        }
    }

    int current_matrix_size = kkt_matrix.matrix.size();
    if (!good_inertia) {
        for (unsigned int i = 0; i < problem.number_variables + this->slacked_constraints.size(); i++) {
            kkt_matrix.add_term(this->inertia_hessian, i, i);
        }
        for (unsigned int j = problem.number_variables + this->slacked_constraints.size(); j < problem.number_variables + this->slacked_constraints.size() + problem.number_constraints; j++) {
            kkt_matrix.add_term(-this->inertia_constraints, j, j);
        }
    }

    while (!good_inertia) {
        std::cout << "Testing factorization with inertia term " << this->inertia_hessian << "\n";
        DEBUG << "KKT matrix:\n";
        for (unsigned int k = 0; k < kkt_matrix.matrix.size(); k++) {
            DEBUG << "m(" << kkt_matrix.row_indices[k] << ", " << kkt_matrix.column_indices[k] << ") = " << kkt_matrix.matrix[k] << "\n";
        }
        this->data = this->solver.factorize(kkt_matrix);

        if (!this->solver.matrix_is_singular() && this->solver.number_negative_eigenvalues() == problem.number_constraints) {
            good_inertia = true;
            std::cout << "Factorization was a success\n";
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
                for (unsigned int i = 0; i < problem.number_variables + this->slacked_constraints.size(); i++) {
                    kkt_matrix.matrix[current_matrix_size + i] = this->inertia_hessian;
                }
                for (unsigned int j = problem.number_variables + this->slacked_constraints.size(); j < problem.number_variables + this->slacked_constraints.size() + problem.number_constraints; j++) {
                    kkt_matrix.matrix[current_matrix_size + j] = -this->inertia_constraints;
                }
            }
        }
    }
    return;
}

std::vector<double> InteriorPoint::generate_rhs(Problem& problem, Iterate& current_iterate, std::vector<double>& variable_lb, std::vector<double>& variable_ub) {
    /* generate the right-hand side */
    std::vector<double> rhs(problem.number_variables + this->slacked_constraints.size() + problem.number_constraints);

    DEBUG << "\n";

    /* objective gradient */
    for (std::pair<int, double> term : current_iterate.objective_gradient) {
        int variable_index = term.first;
        double derivative = term.second;
        DEBUG << "rhs[" << variable_index << "] = " << -problem.objective_sign * derivative << "\n";
        rhs[variable_index] = -problem.objective_sign*derivative;
    }
    /* constraint gradients */
    for (int j = 0; j < problem.number_constraints; j++) {
        double multiplier_j = current_iterate.constraint_multipliers[j];
        if (multiplier_j != 0.) {
            for (std::pair<int, double> term : current_iterate.constraints_jacobian[j]) {
                int variable_index = term.first;
                double derivative = term.second;
                DEBUG << "rhs[" << variable_index << "] += lambda_" << j << "*" << derivative << "\n";
                rhs[variable_index] += multiplier_j*derivative;
            }
        }
    }

    /* variable bound constraints */
    for (int i : this->lower_bounded_variables) {
        DEBUG << "rhs[" << i << "] += " << this->mu << "/(x" << i << " - " << variable_lb[i] << ")\n";
        rhs[i] += this->mu / (current_iterate.x[i] - variable_lb[i]);
    }
    for (int i : this->upper_bounded_variables) {
        DEBUG << "rhs[" << i << "] += " << this->mu << "/(x" << i << " - " << variable_ub[i] << ")\n";
        rhs[i] += this->mu / (current_iterate.x[i] - variable_ub[i]);
    }

    /* slack bound constraints */
    int current_slack = 0;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
            DEBUG << "rhs[" << (problem.number_variables + j) << "] += " << this->mu << "/(s" << current_slack << " - " << problem.constraint_lb[j] << ") - lambda_" << j << "\n";
            rhs[problem.number_variables + current_slack] = this->mu / (current_iterate.x[problem.number_variables + current_slack] - problem.constraint_lb[j]) - current_iterate.constraint_multipliers[j];
        }
        if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
            DEBUG << "rhs[" << (problem.number_variables + j) << "] += " << this->mu << "/(s" << current_slack << " - " << problem.constraint_ub[j] << ") - lambda_" << j << "\n";
            rhs[problem.number_variables + current_slack] = this->mu / (current_iterate.x[problem.number_variables + current_slack] - problem.constraint_ub[j]) - current_iterate.constraint_multipliers[j];
        }
        if (problem.constraint_status[j] != EQUAL_BOUNDS) {
            current_slack++;
        }
    }

    /* constraint evaluations */
    current_slack = 0;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] == EQUAL_BOUNDS) {
            // add the bound
            rhs[problem.number_variables + this->slacked_constraints.size() + j] = -(current_iterate.constraints[j] - problem.constraint_lb[j]);
        }
        else {
            // add the slack
            rhs[problem.number_variables + this->slacked_constraints.size() + j] = -(current_iterate.constraints[j] - current_iterate.x[problem.number_variables + current_slack]);
            current_slack++;
        }
    }
    std::cout << "RHS: ";
    print_vector(std::cout, rhs);

    return rhs;
}

std::vector<double> InteriorPoint::compute_bound_multiplier_displacements(Problem& problem, Iterate& current_iterate, std::vector<double>& solution, std::vector<double>& variable_lb, std::vector<double>& variable_ub) {
    std::vector<double> delta_z(current_iterate.bound_multipliers.size());
    int current_multiplier = 0;
    for (int i : this->lower_bounded_variables) {
        delta_z[current_multiplier] = this->mu / (current_iterate.x[i] - variable_lb[i]) - current_iterate.bound_multipliers[current_multiplier] - current_iterate.bound_multipliers[current_multiplier] / (current_iterate.x[i] - variable_lb[i]) * solution[i];
        current_multiplier++;
    }
    for (int i : this->upper_bounded_variables) {
        delta_z[current_multiplier] = this->mu / (current_iterate.x[i] - variable_ub[i]) - current_iterate.bound_multipliers[current_multiplier] - current_iterate.bound_multipliers[current_multiplier] / (current_iterate.x[i] - variable_ub[i]) * solution[i];
        current_multiplier++;
    }
    int current_slack = 0;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
            delta_z[current_multiplier] = this->mu / (current_iterate.x[problem.number_variables + current_slack] - problem.constraint_lb[j]) - current_iterate.bound_multipliers[current_multiplier] - current_iterate.bound_multipliers[current_multiplier] / (current_iterate.x[problem.number_variables + current_slack] - problem.constraint_lb[j]) * solution[problem.number_variables + current_slack];
            current_multiplier++;
        }
        if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
            delta_z[current_multiplier] = this->mu / (current_iterate.x[problem.number_variables + current_slack] - problem.constraint_ub[j]) - current_iterate.bound_multipliers[current_multiplier] - current_iterate.bound_multipliers[current_multiplier] / (current_iterate.x[problem.number_variables + current_slack] - problem.constraint_ub[j]) * solution[problem.number_variables + current_slack];
            current_multiplier++;
        }
        if (problem.constraint_status[j] != EQUAL_BOUNDS) {
            current_slack++;
        }
    }
    return delta_z;
}

double InteriorPoint::update_barrier_parameter(Problem& problem, Iterate& current_iterate) {
    return this->mu / 10.;
}

void InteriorPoint::compute_measures(Problem& problem, Iterate& iterate) {
    /* evaluate constraints with slacks */
    iterate.feasibility_measure = this->constraint_violation(problem, iterate);
    /* compute barrier objective */
    iterate.optimality_measure = this->objective(problem, iterate);
    return;
}

double InteriorPoint::constraint_violation(Problem& problem, Iterate& iterate) {
    double constraint_violation = 0.;
    int current_slack = 0;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] == EQUAL_BOUNDS) {
            double constraint_value = iterate.constraints[j] - problem.constraint_lb[j];
            constraint_violation += std::abs(constraint_value);
        }
        else {
            double constraint_value = iterate.constraints[j] - iterate.x[problem.number_variables + current_slack];
            constraint_violation += std::abs(constraint_value);
            current_slack++;
        }
    }
    return constraint_violation;
}

double InteriorPoint::objective(Problem& problem, Iterate& iterate) {
    double objective = iterate.objective;
    for (int i : this->lower_bounded_variables) {
        objective -= this->mu * std::log(iterate.x[i] - problem.variable_lb[i]);
    }
    for (int i : this->upper_bounded_variables) {
        objective -= this->mu * std::log(problem.variable_ub[i] - iterate.x[i]);
    }
    int current_slack = 0;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
            objective -= this->mu * std::log(iterate.x[problem.number_variables + current_slack] - problem.constraint_lb[j]);
        }
        if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
            objective -= this->mu * std::log(problem.constraint_ub[j] - iterate.x[problem.number_variables + current_slack]);
        }
        if (problem.constraint_status[j] != EQUAL_BOUNDS) {
            current_slack++;
        }
    }
    return objective;
}

double InteriorPoint::evaluate_local_model(Problem& problem, Iterate& current_iterate, std::vector<double>& solution) {
    double predicted_reduction = 0.;
    current_iterate.compute_objective_gradient(problem);

    /* objective gradient */
    predicted_reduction += dot(solution, current_iterate.objective_gradient);

    /* barrier terms for variables */
    for (int i : this->lower_bounded_variables) {
        predicted_reduction -= this->mu / (current_iterate.x[i] - problem.variable_lb[i]) * solution[i];
    }
    for (int i : this->upper_bounded_variables) {
        predicted_reduction -= this->mu / (current_iterate.x[i] - problem.variable_ub[i]) * solution[i];
    }

    /* barrier terms for slacks */
    int current_slack = 0;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
            predicted_reduction -= this->mu / (current_iterate.x[problem.number_variables + current_slack] - problem.constraint_lb[j]) * solution[problem.number_variables + current_slack];
        }
        if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
            predicted_reduction -= this->mu / (current_iterate.x[problem.number_variables + current_slack] - problem.constraint_ub[j]) * solution[problem.number_variables + current_slack];
        }
        if (problem.constraint_status[j] != EQUAL_BOUNDS) {
            current_slack++;
        }
    }
    return predicted_reduction;
}

double InteriorPoint::compute_error(Problem& problem, Iterate& iterate) {
    double KKT_error = 0.;
    double complementarity_error = 0.;
    double residual = 0.;


    /* residual */
    int current_slack = 0;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] == EQUAL_BOUNDS) {
            residual = std::max(residual, std::abs(iterate.constraints[j] - problem.constraint_lb[j]));
        }
        else {
            residual = std::max(residual, std::abs(iterate.constraints[j] - iterate.x[problem.number_variables + current_slack]));
            current_slack++;
        }
    }

    double error = std::max(KKT_error, std::max(complementarity_error, residual));
    return error;
}

LocalSolution InteriorPoint::compute_infeasibility_step(Problem& problem, Iterate& current_iterate, double radius, LocalSolution& phase_II_solution) {
    std::vector<double> x, z, l;
    return LocalSolution(x, z, l);
}

LocalSolution InteriorPoint::compute_l1_penalty_step(Problem& problem, Iterate& current_iterate, double radius, double penalty_parameter, PenaltyDimensions penalty_dimensions) {
    std::vector<double> x, z, l;
    return LocalSolution(x, z, l);
}