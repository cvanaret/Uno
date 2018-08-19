#include <cmath>
#include "InteriorPoint.hpp"

InteriorPoint::InteriorPoint() : Subproblem("IPM"), mu(0.1), tau_min(0.99), default_multiplier(1.), k_sigma(1e10), inertia_term(0.) {
}

void InteriorPoint::initialize(Problem& problem, Iterate& current_iterate, int number_variables, int number_constraints, bool use_trust_region) {
    /* if trust region is used, bound constraints become range constraints */
    this->variable_status = problem.variable_status;
    if (use_trust_region) {
        for (int i = 0; i < problem.number_variables; i++) {
            if (this->variable_status[i] != EQUAL_BOUNDS) {
                this->variable_status[i] = BOUNDED_BOTH_SIDES;
            }
        }
    }

    /* identify the variable multipliers */
    for (int i = 0; i < problem.number_variables; i++) {
        if (this->variable_status[i] == BOUNDED_LOWER || this->variable_status[i] == BOUNDED_BOTH_SIDES) {
            this->variable_lb.push_back(i);
        }
        if (this->variable_status[i] == BOUNDED_UPPER || this->variable_status[i] == BOUNDED_BOTH_SIDES) {
            this->variable_ub.push_back(i);
        }
    }
    /* identify the inequality constraint slacks */
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] != EQUAL_BOUNDS) {
            this->slacks.push_back(j);
            if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
                this->slack_lb.push_back(j);
            }
            if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
                this->slack_ub.push_back(j);
            }
        }
    }

    /* initialize the slacks and multipliers */
    std::vector<double> bound_multipliers;
    for (int i : this->variable_lb) {
        bound_multipliers.push_back(this->default_multiplier); // positive multiplier
    }
    for (int i : this->variable_ub) {
        bound_multipliers.push_back(-this->default_multiplier); // negative multiplier
    }
    for (int j : this->slack_lb) {
        bound_multipliers.push_back(this->default_multiplier); // positive multiplier
        double slack_value = this->project_variable_in_bounds(current_iterate.constraints[j], problem.constraint_lb[j], problem.constraint_ub[j]);
        current_iterate.x.push_back(slack_value);
    }
    for (int j : this->slack_ub) {
        bound_multipliers.push_back(-this->default_multiplier); // negative multiplier
        double slack_value = this->project_variable_in_bounds(current_iterate.constraints[j], problem.constraint_lb[j], problem.constraint_ub[j]);
        current_iterate.x.push_back(slack_value);
    }
    current_iterate.bound_multipliers = bound_multipliers;

    /* compute first-order information */
    current_iterate.compute_objective_gradient(problem);
    current_iterate.compute_constraints_jacobian(problem);
    /* compute least-square multipliers */
    current_iterate.constraint_multipliers = this->estimate_initial_multipliers(problem, current_iterate);

    std::cout << this->slacks.size() << " slacks\n";
    std::cout << current_iterate.bound_multipliers.size() << " bound multipliers\n";
    std::cout << current_iterate.constraint_multipliers.size() << " constraint multipliers\n";
    std::cout << "variable lb: "; print_vector(std::cout, this->variable_lb);
    std::cout << "variable ub: "; print_vector(std::cout, this->variable_ub);
    std::cout << "slacks: "; print_vector(std::cout, this->slacks);
    std::cout << "slack lb: "; print_vector(std::cout, this->slack_lb);
    std::cout << "slack ub: "; print_vector(std::cout, this->slack_ub);
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


    /* summary */
    for (int i = 0; i < problem.number_variables; i++) {
        std::cout << "x" << i << " in [" << variable_lb[i] << ", " << variable_ub[i] << "]\n";
    }
    std::cout << "\n";
    for (int j = 0; j < problem.number_constraints; j++) {
        std::cout << "c" << j << " in [" << problem.constraint_lb[j] << ", " << problem.constraint_ub[j] << "]\n";
    }
    std::cout << "\n";
    std::cout << "x/s is: ";
    print_vector(std::cout, current_iterate.x);
    std::cout << "λ is: ";
    print_vector(std::cout, current_iterate.constraint_multipliers);
    std::cout << "z is: ";
    print_vector(std::cout, current_iterate.bound_multipliers);
    std::cout << "mu is " << this->mu << "\n";
    std::cout << "Constraints: ";
    print_vector(std::cout, current_iterate.constraints);
    std::cout << "\n";

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

    /* retrieve +Δλ (Nocedal p590) */
    for (int j = 0; j < problem.number_constraints; j++) {
        ipm_solution[problem.number_variables + this->slacks.size() + j] = -ipm_solution[problem.number_variables + this->slacks.size() + j];
    }

    /* compute bound multiplier displacements Δz */
    std::vector<double> delta_z = this->compute_bound_multiplier_displacements(problem, current_iterate, ipm_solution, variable_lb, variable_ub);

    std::cout << "MA57 solution:\n";
    std::cout << "Δx: ";
    print_vector(std::cout, ipm_solution, 0, problem.number_variables);
    std::cout << "Δs: ";
    print_vector(std::cout, ipm_solution, problem.number_variables, this->slacks.size());
    std::cout << "Δz: ";
    print_vector(std::cout, delta_z);
    std::cout << "Δλ: ";
    print_vector(std::cout, ipm_solution, problem.number_variables + this->slacks.size(), problem.number_constraints);

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
    for (unsigned int j = 0; j < this->slacks.size(); j++) {
        trial_x[problem.number_variables + j] = primal_length * ipm_solution[problem.number_variables + j];
    }
    std::cout << "Scaled Δx/Δs: ";
    print_vector(std::cout, trial_x);

    for (unsigned int j = 0; j < current_iterate.bound_multipliers.size(); j++) {
        trial_bound_multipliers[j] = current_iterate.bound_multipliers[j] + dual_length * delta_z[j];
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
        trial_constraint_multipliers[j] = primal_length * ipm_solution[problem.number_variables + this->slacks.size() + j];
    }
    std::cout << "Scaled Δλ: ";
    print_vector(std::cout, trial_constraint_multipliers);

    LocalSolution solution(trial_x, trial_bound_multipliers, trial_constraint_multipliers);
    // TODO
    solution.norm = 1234.;
    return solution;
}

std::vector<double> InteriorPoint::estimate_initial_multipliers(Problem& problem, Iterate& current_iterate) {
    // TODO handle number of multipliers

    /* build the symmetric matrix */
    COOMatrix matrix(problem.number_variables + problem.number_constraints, 0);

    /* identity blocks */
    for (int i = 0; i < problem.number_variables; i++) {
        matrix.add_term(-1., i, i);
    }
    for (unsigned int current_slack = 0; current_slack < this->slacks.size(); current_slack++) {
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
    for (int i : this->variable_lb) {
        rhs[i] -= current_iterate.bound_multipliers[current_multiplier];
        current_multiplier++;
    }
    for (int i : this->variable_ub) {
        rhs[i] -= current_iterate.bound_multipliers[current_multiplier];
        current_multiplier++;
    }

    /* slack bound constraints */
    for (int j : this->slack_lb) {
        DEBUG << "rhs[" << (problem.number_variables + j) << "] += " << current_iterate.bound_multipliers[current_multiplier] << "\n";
        rhs[problem.number_variables + j] += current_iterate.bound_multipliers[current_multiplier];
        current_multiplier++;
    }
    for (int j : this->slack_ub) {
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
    for (int i : this->variable_lb) {
        double trial_alpha_xi = -tau * (current_iterate.x[i] - variable_lb[i]) / ipm_solution[i];
        if (0 < trial_alpha_xi && trial_alpha_xi <= 1.) {
            primal_length = std::min(primal_length, trial_alpha_xi);
        }
    }
    for (int i : this->variable_ub) {
        double trial_alpha_xi = -tau * (current_iterate.x[i] - variable_ub[i]) / ipm_solution[i];
        if (0 < trial_alpha_xi && trial_alpha_xi <= 1.) {
            primal_length = std::min(primal_length, trial_alpha_xi);
        }
    }
    int current_slack = 0;
    for (int j : this->slack_lb) {
        double trial_primal_lengthj = -tau * (current_iterate.x[problem.number_variables + current_slack] - problem.constraint_lb[j]) / ipm_solution[problem.number_variables + current_slack];
        if (0 < trial_primal_lengthj && trial_primal_lengthj <= 1.) {
            primal_length = std::min(primal_length, trial_primal_lengthj);
        }
        current_slack++;
    }
    for (int j : this->slack_ub) {
        double trial_primal_lengthj = -tau * (current_iterate.x[problem.number_variables + current_slack] - problem.constraint_ub[j]) / ipm_solution[problem.number_variables + current_slack];
        if (0 < trial_primal_lengthj && trial_primal_lengthj <= 1.) {
            primal_length = std::min(primal_length, trial_primal_lengthj);
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

LocalSolution InteriorPoint::compute_infeasibility_step(Problem& problem, Iterate& current_iterate, double radius, LocalSolution& phase_II_solution) {
    std::vector<double> x, z, l;
    return LocalSolution(x, z, l);
}

LocalSolution InteriorPoint::compute_l1_penalty_step(Problem& problem, Iterate& current_iterate, double radius, double penalty_parameter, PenaltyDimensions penalty_dimensions) {
    std::vector<double> x, z, l;
    return LocalSolution(x, z, l);
}

COOMatrix InteriorPoint::generate_kkt_matrix(Problem& problem, Iterate& current_iterate, std::vector<double>& variable_lb, std::vector<double>& variable_ub) {
    /* compute the Lagrangian Hessian */
    current_iterate.compute_hessian(problem, problem.objective_sign, current_iterate.constraint_multipliers);
    COOMatrix kkt_matrix = current_iterate.hessian.to_COO();
    kkt_matrix.size = problem.number_variables + this->slacks.size() + problem.number_constraints;

    /* variable bound constraints */
    int current_multiplier = 0;
    for (int i : this->variable_lb) {
        kkt_matrix.add_term(current_iterate.bound_multipliers[current_multiplier] / (current_iterate.x[i] - variable_lb[i]), i, i);
        current_multiplier++;
    }
    for (int i : this->variable_ub) {
        kkt_matrix.add_term(current_iterate.bound_multipliers[current_multiplier] / (current_iterate.x[i] - variable_ub[i]), i, i);
        current_multiplier++;
    }

    /* slack bound constraints */
    int current_column = problem.number_variables;
    int current_slack = 0;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
            kkt_matrix.add_term(current_iterate.bound_multipliers[current_multiplier] / (current_iterate.x[problem.number_variables + current_slack] - problem.constraint_lb[j]), current_column, current_column);
            kkt_matrix.add_term(-1., current_column, problem.number_variables + this->slacks.size() + j);
            current_multiplier++;
        }
        if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
            kkt_matrix.add_term(current_iterate.bound_multipliers[current_multiplier] / (current_iterate.x[problem.number_variables + current_slack] - problem.constraint_ub[j]), current_column, current_column);
            kkt_matrix.add_term(-1., current_column, problem.number_variables + this->slacks.size() + j);
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
    int current_matrix_size = kkt_matrix.matrix.size();
    for (unsigned int i = 0; i < problem.number_variables + this->slacks.size(); i++) {
        kkt_matrix.add_term(this->inertia_term, i, i);
    }

    bool good_inertia = false;
    while (!good_inertia) {
        std::cout << "Testing factorization with inertia term " << this->inertia_term << "\n";
        DEBUG << "KKT matrix:\n";
        for (unsigned int k = 0; k < kkt_matrix.matrix.size(); k++) {
            DEBUG << "m(" << kkt_matrix.row_indices[k] << ", " << kkt_matrix.column_indices[k] << ") = " << kkt_matrix.matrix[k] << "\n";
        }
        this->data = this->solver.factorize(kkt_matrix);
        std::cout << "The matrix has " << this->solver.number_negative_eigenvalues() << " negative eigenvalues\n";

        if (this->solver.number_negative_eigenvalues() == problem.number_constraints) {
            good_inertia = true;
            std::cout << "Factorization was a success\n";
        }
        else {
            std::cout << "Bad inertia with inertia term " << this->inertia_term << "\n";
            /* increase inertia factor */
            this->inertia_term = (this->inertia_term == 0.) ? 1e-4 : 100 * this->inertia_term;
            for (unsigned int i = 0; i < problem.number_variables + this->slacks.size(); i++) {
                kkt_matrix.matrix[current_matrix_size + i] = this->inertia_term;
            }
        }
    }

    return kkt_matrix;
}

std::vector<double> InteriorPoint::generate_rhs(Problem& problem, Iterate& current_iterate, std::vector<double>& variable_lb, std::vector<double>& variable_ub) {
    /* generate the right-hand side */
    std::vector<double> rhs(problem.number_variables + this->slacks.size() + problem.number_constraints);

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
                DEBUG << "rhs[" << variable_index << "] += " << multiplier_j * derivative << "\n";
                rhs[variable_index] += multiplier_j*derivative;
            }
        }
    }

    /* variable bound constraints */
    for (int i : this->variable_lb) {
        DEBUG << "rhs[" << i << "] += " << this->mu << "/(x" << i << " - " << variable_lb[i] << ")\n";
        rhs[i] += this->mu / (current_iterate.x[i] - variable_lb[i]);
    }
    for (int i : this->variable_ub) {
        DEBUG << "rhs[" << i << "] += " << this->mu << "/(x" << i << " - " << variable_ub[i] << ")\n";
        rhs[i] += this->mu / (current_iterate.x[i] - variable_ub[i]);
    }

    /* slack bound constraints */
    int current_slack = 0;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
            DEBUG << "rhs[" << (problem.number_variables + j) << "] += " << this->mu << "/(s" << current_slack << " - " << problem.constraint_lb[j] << ") - " << current_iterate.constraint_multipliers[j] << "\n";
            rhs[problem.number_variables + current_slack] = this->mu / (current_iterate.x[problem.number_variables + current_slack] - problem.constraint_lb[j]) - current_iterate.constraint_multipliers[j];
        }
        if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
            DEBUG << "rhs[" << (problem.number_variables + j) << "] += " << this->mu << "/(s" << current_slack << " - " << problem.constraint_ub[j] << ") - " << current_iterate.constraint_multipliers[j] << "\n";
            rhs[problem.number_variables + current_slack] = this->mu / (current_iterate.x[problem.number_variables + current_slack] - problem.constraint_ub[j]) - current_iterate.constraint_multipliers[j];
        }
        if (problem.constraint_status[j] != EQUAL_BOUNDS) {
            current_slack++;
        }
    }

    /* constraint evaluations */
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] == EQUAL_BOUNDS) {
            rhs[problem.number_variables + this->slacks.size() + j] = -(current_iterate.constraints[j] - problem.constraint_lb[j]);
        }
        else {
            rhs[problem.number_variables + this->slacks.size() + j] = -current_iterate.constraints[j];
        }
    }
    /* if inequality constraint, add corresponding slack */
    for (unsigned int current_slack = 0; current_slack < this->slacks.size(); current_slack++) {
        int j = this->slacks[current_slack];
        rhs[problem.number_variables + this->slacks.size() + j] += current_iterate.x[problem.number_variables + current_slack];
        DEBUG << "rhs[" << (problem.number_variables + this->slacks.size() + j) << "] = " << rhs[problem.number_variables + this->slacks.size() + j] << "\n";
    }
    std::cout << "RHS: ";
    print_vector(std::cout, rhs);

    return rhs;
}

std::vector<double> InteriorPoint::compute_bound_multiplier_displacements(Problem& problem, Iterate& current_iterate, std::vector<double>& solution, std::vector<double>& variable_lb, std::vector<double>& variable_ub) {
    std::vector<double> delta_z(current_iterate.bound_multipliers.size());
    int current_multiplier = 0;
    for (int i : this->variable_lb) {
        delta_z[current_multiplier] = this->mu / (current_iterate.x[i] - variable_lb[i]) - current_iterate.bound_multipliers[current_multiplier] - current_iterate.bound_multipliers[current_multiplier] / (current_iterate.x[i] - variable_lb[i]) * solution[i];
        current_multiplier++;
    }
    for (int i : this->variable_ub) {
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
