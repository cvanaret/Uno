#include <cmath>
#include "InteriorPoint.hpp"
#include "Argonot.hpp"

InteriorPoint::InteriorPoint(Problem& problem, std::string hessian_evaluation_method, bool use_trust_region, bool scale_residuals):
Subproblem("l2", problem.variables_bounds, scale_residuals), // use the l2 norm to compute residuals
hessian_evaluation(HessianEvaluationFactory::create(hessian_evaluation_method, problem.number_variables)),
mu_optimality(0.1), mu_feasibility(mu_optimality), inertia_hessian(0.), inertia_hessian_last(0.),
inertia_constraints(0.), default_multiplier(1.), iteration(0),
parameters({0.99, 1e10, 100., 0.2, 1.5, 10., 1e10}) {
    /* the subproblem optimizes the original variables and slacks for inequality constraints */
    int number_variables = problem.number_variables + problem.inequality_constraints.size();
    this->subproblem_variables_bounds.resize(number_variables);
    
    /* identify the bounded variables */
    for (int i = 0; i < problem.number_variables; i++) {
        if (use_trust_region || (problem.variable_status[i] == BOUNDED_LOWER || problem.variable_status[i] == BOUNDED_BOTH_SIDES)) {
            this->lower_bounded_variables.insert(i);
        }
        if (use_trust_region || (problem.variable_status[i] == BOUNDED_UPPER || problem.variable_status[i] == BOUNDED_BOTH_SIDES)) {
            this->upper_bounded_variables.insert(i);
        }
    }
    /* identify the inequality constraint slacks */
    for (std::pair<const int, int>& element: problem.inequality_constraints) {
        int j = element.first;
        int slack_index = problem.number_variables + element.second;
        if (problem.constraint_status[j] == BOUNDED_LOWER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
            this->lower_bounded_variables.insert(slack_index);
        }
        if (problem.constraint_status[j] == BOUNDED_UPPER || problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
            this->upper_bounded_variables.insert(slack_index);
        }
        // register the bounds of the slacks
        this->subproblem_variables_bounds[slack_index] = problem.constraint_bounds[j];
    }
}

Iterate InteriorPoint::evaluate_initial_point(Problem& problem, std::vector<double>& x, Multipliers& default_multipliers) {
    int number_variables = problem.number_variables + problem.inequality_constraints.size();
    
    /* make the initial point strictly feasible */
    std::vector<double> reformulated_x(number_variables);
    for (int i = 0; i < problem.number_variables; i++) {
        reformulated_x[i] = Subproblem::project_variable_in_bounds(x[i], problem.variables_bounds[i]);
    }
    
    Multipliers multipliers(number_variables, problem.number_constraints);
    /* generate the bound multipliers */
    for (int i: this->lower_bounded_variables) {
        multipliers.lower_bounds[i] = this->default_multiplier; // positive multiplier
    }
    for (int i: this->upper_bounded_variables) {
        multipliers.upper_bounds[i] = -this->default_multiplier; // negative multiplier
    }
    
    /* generate the first iterate */
    Iterate first_iterate(reformulated_x, multipliers);

    /* initialize the slacks */
    first_iterate.compute_constraints(problem);
    for (std::pair<const int, int>& element: problem.inequality_constraints) {
        int j = element.first;
        int slack_index = problem.number_variables + element.second;
        double slack_value = Subproblem::project_variable_in_bounds(first_iterate.constraints[j], problem.constraint_bounds[j]);
        first_iterate.x[slack_index] = slack_value;
    }

    /* evaluate the constraint Jacobian */
    first_iterate.compute_constraints_jacobian(problem);
    // contribution of the slacks
    for (std::pair<const int, int>& element: problem.inequality_constraints) {
        int j = element.first;
        int slack_index = problem.number_variables + element.second;
        first_iterate.constraints_jacobian[j][slack_index] = -1.;
    }
    /* compute least-square multipliers */
    if (0 < problem.number_constraints) {
        first_iterate.multipliers.constraints = Subproblem::compute_least_square_multipliers(problem, first_iterate, default_multipliers.constraints, this->solver);
    }

    DEBUG << problem.inequality_constraints.size() << " slacks\n";
    DEBUG << first_iterate.multipliers.lower_bounds.size() << " bound multipliers\n";
    DEBUG << first_iterate.multipliers.constraints.size() << " constraint multipliers\n";
    DEBUG << "variable lb: "; print_vector(DEBUG, this->lower_bounded_variables);
    DEBUG << "variable ub: "; print_vector(DEBUG, this->upper_bounded_variables);

    /* compute the optimality and feasibility measures of the initial point */
    this->compute_optimality_measures(problem, first_iterate);
    
    /* if no trust region is used, the problem should be convexified. The inertia of the augmented matrix will be corrected later */
    this->hessian_evaluation->convexify = false;
    
    return first_iterate;
}

double InteriorPoint::compute_KKT_error_scaling(Iterate& current_iterate) {
    /* KKT error */
    double norm_1_constraint_multipliers = norm_1(current_iterate.multipliers.constraints);
    double norm_1_bound_multipliers = norm_1(current_iterate.multipliers.lower_bounds) + norm_1(current_iterate.multipliers.upper_bounds);
    double sd = std::max(this->parameters.smax, (norm_1_constraint_multipliers + norm_1_bound_multipliers) / (current_iterate.x.size() + current_iterate.multipliers.constraints.size())) / this->parameters.smax;
    return sd;
}

/* reduced primal-dual approach */
SubproblemSolution InteriorPoint::compute_optimality_step(Problem& problem, Iterate& current_iterate, double /*trust_region_radius*/) {
    DEBUG << "\nCurrent iterate: ";
    DEBUG << current_iterate;

    current_iterate.compute_constraints_jacobian(problem);
    // contribution of the slacks
    for (std::pair<const int, int>& element: problem.inequality_constraints) {
        int j = element.first;
        int slack_index = problem.number_variables + element.second;
        current_iterate.constraints_jacobian[j][slack_index] = -1.;
    }
    
    /* scaled error terms */
    this->compute_residuals(problem, current_iterate, current_iterate.multipliers, 1.);
    double sd = this->compute_KKT_error_scaling(current_iterate);
    double KKTerror = current_iterate.residuals.KKT / sd;
    double central_complementarity_error = this->compute_central_complementarity_error(current_iterate, this->mu_optimality, this->subproblem_variables_bounds);
    DEBUG << "IPM error (KKT: " << KKTerror << ", cmpl: " << central_complementarity_error << ", feas: " << current_iterate.residuals.constraints << ")\n";

    /* update of the barrier problem */
    double error = std::max(KKTerror, std::max(central_complementarity_error, current_iterate.residuals.constraints));
    if (error <= this->parameters.k_epsilon * this->mu_optimality) {
        // TODO pass tolerance
        double tolerance = 1e-8;
        this->mu_optimality = std::max(tolerance / 10., std::min(this->parameters.k_mu * this->mu_optimality, std::pow(this->mu_optimality, this->parameters.theta_mu)));
        DEBUG << "IPM: mu updated to " << this->mu_optimality << " and filter reset\n";
        // signal the redefinition of the problem to the globalization strategy
        this->subproblem_definition_changed = true;
    }
    DEBUG << "mu is " << this->mu_optimality << "\n";
    this->iteration++;

    /* evaluate the functions at the current iterate */
    this->evaluate_optimality_iterate(problem, current_iterate);
    
    /************************/
    /* solve the KKT system */
    /************************/
    /* KKT matrix */
    COOMatrix kkt_matrix = this->generate_optimality_kkt_matrix(problem, current_iterate, this->subproblem_variables_bounds);

    /* inertia correction */
    MA57Factorization factorization = this->modify_inertia(kkt_matrix, current_iterate.x.size(), problem.number_constraints);
    DEBUG << "KKT matrix:\n" << kkt_matrix << "\n";

    /* right-hand side */
    std::vector<double> rhs = this->generate_kkt_rhs(problem, current_iterate);

    /* compute the solution (Δx, -Δλ) */
    this->solver.solve(factorization, rhs);
    this->number_subproblems_solved++;
    std::vector<double>& solution_IPM = rhs;

    /* generate IPM direction */
    SubproblemSolution solution = this->generate_direction(problem, current_iterate, solution_IPM);
    solution.status = OPTIMAL;
    solution.norm = norm_inf(solution.x, problem.number_variables);
    solution.predicted_reduction = [&](double step_length) {
        return this->compute_predicted_reduction(solution, step_length);
    };

    /* evaluate the barrier objective */
    solution.objective = this->evaluate_local_model(problem, current_iterate, solution.x);
    return solution;
}

void InteriorPoint::evaluate_optimality_iterate(Problem& problem, Iterate& current_iterate) {
    /* compute barrier gradient */
    current_iterate.compute_objective_gradient(problem);
    // contribution of bound constraints
    for (int i: this->lower_bounded_variables) {
        current_iterate.objective_gradient[i] -= this->mu_optimality / (current_iterate.x[i] - this->subproblem_variables_bounds[i].lb);
    }
    for (int i: this->upper_bounded_variables) {
        current_iterate.objective_gradient[i] -= this->mu_optimality / (current_iterate.x[i] - this->subproblem_variables_bounds[i].ub);
    }

    /* compute constraint Jacobian */
    current_iterate.compute_constraints_jacobian(problem);
    // contribution of the slacks
    for (std::pair<const int, int>& element: problem.inequality_constraints) {
        int j = element.first;
        int slack_index = problem.number_variables + element.second;
        current_iterate.constraints_jacobian[j][slack_index] = -1.;
    }
    
    /* compute second-order information */
    this->hessian_evaluation->compute(problem, current_iterate, problem.objective_sign, current_iterate.multipliers.constraints);
    return;
}

SubproblemSolution InteriorPoint::generate_direction(Problem& problem, Iterate& current_iterate, std::vector<double>& solution_IPM) {
    int number_variables = problem.number_variables + problem.inequality_constraints.size();

    /* retrieve +Δλ (Nocedal p590) */
    for (int j = 0; j < problem.number_constraints; j++) {
        int multiplier_index = number_variables + j;
        solution_IPM[multiplier_index] = -solution_IPM[multiplier_index];
    }

    /* compute bound multiplier displacements Δz */
    std::vector<double> lower_delta_z = this->compute_lower_bound_multiplier_displacements(current_iterate, solution_IPM, subproblem_variables_bounds, this->mu_optimality);
    std::vector<double> upper_delta_z = this->compute_upper_bound_multiplier_displacements(current_iterate, solution_IPM, subproblem_variables_bounds, this->mu_optimality);

    /* "fraction to boundary" rule for variables and bound multipliers */
    std::vector<double> trial_x(current_iterate.x.size());
    Multipliers trial_multipliers(current_iterate.x.size(), problem.number_constraints);
    double tau = std::max(this->parameters.tau_min, 1. - this->mu_optimality);
    // scale primal variables and constraints multipliers
    double primal_length = this->compute_primal_length(current_iterate, solution_IPM, subproblem_variables_bounds, tau);
    for (int i = 0; i < number_variables; i++) {
        trial_x[i] = primal_length * solution_IPM[i];
    }
    for (int j = 0; j < problem.number_constraints; j++) {
        trial_multipliers.constraints[j] = current_iterate.multipliers.constraints[j] + primal_length * solution_IPM[number_variables + j];
    }

    // scale dual variables
    double dual_length = this->compute_dual_length(current_iterate, tau, lower_delta_z, upper_delta_z);
    for (int i = 0; i < number_variables; i++) {
        trial_multipliers.lower_bounds[i] = current_iterate.multipliers.lower_bounds[i] + dual_length * lower_delta_z[i];
        trial_multipliers.upper_bounds[i] = current_iterate.multipliers.upper_bounds[i] + dual_length * upper_delta_z[i];
        // rescale the multipliers (IPOPT paper p6)
        //trial_multipliers.lower_bounds[i] = std::max(std::min(trial_multipliers.lower_bounds[i], this->parameters.kappa*this->mu_optimality/(trial_x[i] - subproblem_variables_bounds[i].lb)), this->mu_optimality/(this->parameters.kappa*(trial_x[i] - subproblem_variables_bounds[i].lb)));
    }

    DEBUG << "MA57 solution:\n";
    DEBUG << "Δx: "; print_vector(DEBUG, solution_IPM, 0, problem.number_variables);
    DEBUG << "Δs: "; print_vector(DEBUG, solution_IPM, problem.number_variables, problem.inequality_constraints.size());
    DEBUG << "Δz_L: "; print_vector(DEBUG, lower_delta_z);
    DEBUG << "Δz_U: "; print_vector(DEBUG, upper_delta_z);
    DEBUG << "Δλ: "; print_vector(DEBUG, solution_IPM, number_variables, problem.number_constraints);
    DEBUG << "primal length = " << primal_length << "\n";
    DEBUG << "dual length = " << dual_length << "\n\n";

    SubproblemSolution solution(trial_x, trial_multipliers);
    return solution;
}

double InteriorPoint::compute_primal_length(Iterate& current_iterate, std::vector<double>& ipm_solution, std::vector<Range>& variables_bounds, double tau) {
    double primal_length = 1.;
    for (int i: this->lower_bounded_variables) {
        double trial_alpha_xi = -tau * (current_iterate.x[i] - variables_bounds[i].lb) / ipm_solution[i];
        if (0 < trial_alpha_xi && trial_alpha_xi <= 1.) {
            primal_length = std::min(primal_length, trial_alpha_xi);
        }
    }
    for (int i: this->upper_bounded_variables) {
        double trial_alpha_xi = -tau * (current_iterate.x[i] - variables_bounds[i].ub) / ipm_solution[i];
        if (0 < trial_alpha_xi && trial_alpha_xi <= 1.) {
            primal_length = std::min(primal_length, trial_alpha_xi);
        }
    }
    return primal_length;
}

double InteriorPoint::compute_dual_length(Iterate& current_iterate, double tau, std::vector<double>& lower_delta_z, std::vector<double>& upper_delta_z) {
    double dual_length = 1.;
    for (unsigned int i = 0; i < current_iterate.multipliers.lower_bounds.size(); i++) {
        double trial_alpha_zj = -tau * current_iterate.multipliers.lower_bounds[i] / lower_delta_z[i];
        if (0 < trial_alpha_zj && trial_alpha_zj <= 1.) {
            dual_length = std::min(dual_length, trial_alpha_zj);
        }
        trial_alpha_zj = -tau * current_iterate.multipliers.upper_bounds[i] / upper_delta_z[i];
        if (0 < trial_alpha_zj && trial_alpha_zj <= 1.) {
            dual_length = std::min(dual_length, trial_alpha_zj);
        }
    }
    return dual_length;
}

COOMatrix InteriorPoint::generate_optimality_kkt_matrix(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds) {
    int number_variables = problem.number_variables + problem.inequality_constraints.size();

    /* compute the Lagrangian Hessian */
    COOMatrix kkt_matrix = current_iterate.hessian.to_COO();
    kkt_matrix.dimension = number_variables + problem.number_constraints;

    /* variable bound constraints */
    for (int i: this->lower_bounded_variables) {
        kkt_matrix.add_term(current_iterate.multipliers.lower_bounds[i] / (current_iterate.x[i] - variables_bounds[i].lb), i, i);
    }
    for (int i: this->upper_bounded_variables) {
        kkt_matrix.add_term(current_iterate.multipliers.upper_bounds[i] / (current_iterate.x[i] - variables_bounds[i].ub), i, i);
    }

    /* Jacobian of general constraints */
    for (int j = 0; j < problem.number_constraints; j++) {
        for (std::pair<int, double> term: current_iterate.constraints_jacobian[j]) {
            int variable_index = term.first;
            double derivative = term.second;
            kkt_matrix.add_term(derivative, variable_index, number_variables + j);
        }
    }
    return kkt_matrix;
}

MA57Factorization InteriorPoint::modify_inertia(COOMatrix& kkt_matrix, int size_first_block, int size_second_block) {
    this->inertia_hessian = 0.;
    this->inertia_constraints = 0.;
    DEBUG << "Testing factorization with inertia term " << this->inertia_hessian << "\n";
    MA57Factorization factorization = this->solver.factorize(kkt_matrix);

    bool good_inertia = false;
    if (!factorization.matrix_is_singular() && factorization.number_negative_eigenvalues() == size_second_block) {
        DEBUG << "Factorization was a success\n";
        good_inertia = true;
    }
    else {
        // inertia term for constraints
        if (factorization.matrix_is_singular()) {
            DEBUG << "Matrix is singular\n";
            this->inertia_constraints = 1e-8 * std::pow(this->mu_optimality, 0.25);
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
        for (int i = 0; i < size_first_block; i++) {
            kkt_matrix.add_term(this->inertia_hessian, i, i);
        }
        for (int j = size_first_block; j < size_first_block + size_second_block; j++) {
            kkt_matrix.add_term(-this->inertia_constraints, j, j);
        }
    }

    while (!good_inertia) {
        DEBUG << "Testing factorization with inertia term " << this->inertia_hessian << "\n";
        factorization = this->solver.factorize(kkt_matrix);

        if (!factorization.matrix_is_singular() && factorization.number_negative_eigenvalues() == size_second_block) {
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
                throw UnstableInertiaCorrection();
            }
            else {
                for (int i = 0; i < size_first_block; i++) {
                    kkt_matrix.matrix[current_matrix_size + i] = this->inertia_hessian;
                }
                for (int j = size_first_block; j < size_first_block + size_second_block; j++) {
                    kkt_matrix.matrix[current_matrix_size + j] = -this->inertia_constraints;
                }
            }
        }
    }
    return factorization;
}

std::vector<double> InteriorPoint::generate_kkt_rhs(Problem& problem, Iterate& current_iterate) {
    int number_variables = problem.number_variables + problem.inequality_constraints.size();

    /* generate the right-hand side */
    std::vector<double> rhs(number_variables + problem.number_constraints);

    /* barrier objective gradient */
    for (std::pair<int, double> term: current_iterate.objective_gradient) {
        int i = term.first;
        double derivative = term.second;
        rhs[i] = -derivative;
    }

    /* constraint gradients */
    for (int j = 0; j < problem.number_constraints; j++) {
        if (current_iterate.multipliers.constraints[j] != 0.) {
            for (std::pair<int, double> term: current_iterate.constraints_jacobian[j]) {
                int variable_index = term.first;
                double derivative = term.second;
                rhs[variable_index] += current_iterate.multipliers.constraints[j] * derivative;
            }
        }
    }

    /* constraint evaluations */
    int slack_index = 0;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] == EQUAL_BOUNDS) {
            // add the bound
            rhs[number_variables + j] = -(current_iterate.constraints[j] - problem.constraint_bounds[j].lb);
        }
        else {
            // add the slack
            rhs[number_variables + j] = -(current_iterate.constraints[j] - current_iterate.x[problem.number_variables + slack_index]);
            slack_index++;
        }
    }
    DEBUG << "RHS: ";
    print_vector(DEBUG, rhs);

    return rhs;
}

std::vector<double> InteriorPoint::compute_lower_bound_multiplier_displacements(Iterate& current_iterate, std::vector<double>& solution, std::vector<Range>& variables_bounds, double mu) {
    std::vector<double> delta_z(current_iterate.multipliers.lower_bounds.size());
    for (int i: this->lower_bounded_variables) {
        delta_z[i] = mu / (current_iterate.x[i] - variables_bounds[i].lb) - current_iterate.multipliers.lower_bounds[i] - current_iterate.multipliers.lower_bounds[i] / (current_iterate.x[i] - variables_bounds[i].lb) * solution[i];
    }
    return delta_z;
}

std::vector<double> InteriorPoint::compute_upper_bound_multiplier_displacements(Iterate& current_iterate, std::vector<double>& solution, std::vector<Range>& variables_bounds, double mu) {
    std::vector<double> delta_z(current_iterate.multipliers.upper_bounds.size());
    for (int i: this->upper_bounded_variables) {
        delta_z[i] = mu / (current_iterate.x[i] - variables_bounds[i].ub) - current_iterate.multipliers.upper_bounds[i] - current_iterate.multipliers.upper_bounds[i] / (current_iterate.x[i] - variables_bounds[i].ub) * solution[i];
    }
    return delta_z;
}

void InteriorPoint::compute_optimality_measures(Problem& problem, Iterate& iterate) {
    /* evaluate constraints with slacks */
    iterate.feasibility_measure = this->constraint_violation(problem, iterate);
    /* compute barrier objective */
    iterate.optimality_measure = this->barrier_function(problem, iterate, this->subproblem_variables_bounds);
    return;
}

void InteriorPoint::compute_infeasibility_measures(Problem& problem, Iterate& iterate, SubproblemSolution& /*solution*/) {
    this->compute_optimality_measures(problem, iterate);
    return;
}

double InteriorPoint::constraint_violation(Problem& problem, Iterate& iterate) {
    iterate.compute_constraints(problem);
    // compute l2 square norm
    std::vector<double> residuals(problem.number_constraints);
    int slack_index = problem.number_variables;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] == EQUAL_BOUNDS) {
            double constraint_value = iterate.constraints[j] - problem.constraint_bounds[j].lb;
            residuals[j] = constraint_value;
        }
        else {
            double constraint_value = iterate.constraints[j] - iterate.x[slack_index];
            residuals[j] = constraint_value;
            slack_index++;
        }
    }
    return norm(residuals, this->residual_norm);
}

double InteriorPoint::barrier_function(Problem& problem, Iterate& iterate, std::vector<Range>& variables_bounds) {
    /* original objective */
    iterate.compute_objective(problem);
    double objective = iterate.objective;

    /* bound constraints */
    for (int i: this->lower_bounded_variables) {
        objective -= this->mu_optimality * std::log(iterate.x[i] - variables_bounds[i].lb);
    }
    for (int i: this->upper_bounded_variables) {
        objective -= this->mu_optimality * std::log(variables_bounds[i].ub - iterate.x[i]);
    }
    return objective;
}

double InteriorPoint::evaluate_local_model(Problem& /*problem*/, Iterate& current_iterate, std::vector<double>& solution) {
    double subproblem_objective = dot(solution, current_iterate.objective_gradient);
    return subproblem_objective;
}

double InteriorPoint::compute_predicted_reduction(SubproblemSolution& solution, double step_length) {
    // the predicted reduction is linear
    return -step_length*solution.objective;
}

SubproblemSolution InteriorPoint::compute_infeasibility_step(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution, double trust_region_radius) {
    int number_variables = problem.number_variables + problem.inequality_constraints.size();

    DEBUG << "restoration x: "; print_vector(DEBUG, current_iterate.x);
    
    /* multipliers = 2*c */
    std::vector<double> restoration_multipliers(problem.number_constraints);
    int slack_index = 0;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] == EQUAL_BOUNDS) {
            // add the bound
            restoration_multipliers[j] = 2*(current_iterate.constraints[j] - problem.constraint_bounds[j].lb);
        }
        else {
            // add the slack
            restoration_multipliers[j] = 2*(current_iterate.constraints[j] - current_iterate.x[problem.number_variables + slack_index]);
            slack_index++;
        }
    }
    
    /* compute the Lagrangian Hessian */
    this->hessian_evaluation->compute(problem, current_iterate, 0., restoration_multipliers);
    ArgonotMatrix kkt_matrix = current_iterate.hessian.to_ArgonotMatrix(number_variables);
    // contribution of 2 \nabla c \nabla c^T
    for (int j = 0; j < problem.number_constraints; j++) {
        DEBUG << "Gradient c" << j << ": "; print_vector(DEBUG, current_iterate.constraints_jacobian[j]);
        kkt_matrix.add_outer_product(current_iterate.constraints_jacobian[j], 2.);
    }
    // variable bound constraints
    for (int i: this->lower_bounded_variables) {
        kkt_matrix.add_term(current_iterate.multipliers.lower_bounds[i] / (current_iterate.x[i] - this->subproblem_variables_bounds[i].lb), i, i);
    }
    for (int i: this->upper_bounded_variables) {
        kkt_matrix.add_term(current_iterate.multipliers.upper_bounds[i] / (current_iterate.x[i] - this->subproblem_variables_bounds[i].ub), i, i);
    }
    
    /* factorization by the linear solver */
    COOMatrix coo_matrix = kkt_matrix.to_COO();
    
    /* inertia correction */
    MA57Factorization factorization = this->modify_inertia(coo_matrix, current_iterate.x.size(), 0);
    
    DEBUG << "restoration KKT matrix:\n" << coo_matrix;
    
    /* right-hand side */
    std::vector<double> rhs(number_variables);
    // constraint Jacobian
    for (int j = 0; j < problem.number_constraints; j++) {
        if (restoration_multipliers[j] != 0.) {
            for (std::pair<int, double> term: current_iterate.constraints_jacobian[j]) {
                int i = term.first;
                double derivative = term.second;
                rhs[i] += restoration_multipliers[j] * derivative;
            }
        }
    }
    // variable bound constraints
    for (int i: this->lower_bounded_variables) {
        rhs[i] += this->mu_feasibility/ (current_iterate.x[i] - this->subproblem_variables_bounds[i].lb);
    }
    for (int i: this->upper_bounded_variables) {
        rhs[i] += this->mu_feasibility / (current_iterate.x[i] - this->subproblem_variables_bounds[i].ub);
    }
    DEBUG << "restoration RHS: "; print_vector(DEBUG, rhs); DEBUG << "\n";
    
    /* compute the solution Δx */
    this->solver.solve(factorization, rhs);
    this->number_subproblems_solved++;
    std::vector<double>& solution_IPM = rhs;
    
    /* compute bound multiplier displacements Δz */
    std::vector<double> lower_delta_z = this->compute_lower_bound_multiplier_displacements(current_iterate, solution_IPM, subproblem_variables_bounds, this->mu_feasibility);
    std::vector<double> upper_delta_z = this->compute_upper_bound_multiplier_displacements(current_iterate, solution_IPM, subproblem_variables_bounds, this->mu_feasibility);
    
    /* create the solution */
    std::vector<double> trial_x(current_iterate.x.size());
    Multipliers trial_multipliers(current_iterate.x.size(), current_iterate.constraints.size());
    double tau = std::max(this->parameters.tau_min, 1. - this->mu_feasibility);
    // scale primal variables and constraints multipliers
    double primal_length = this->compute_primal_length(current_iterate, solution_IPM, subproblem_variables_bounds, tau);
    for (int i = 0; i < number_variables; i++) {
        trial_x[i] = primal_length * solution_IPM[i];
    }
    // scale dual variables
    double dual_length = this->compute_dual_length(current_iterate, tau, lower_delta_z, upper_delta_z);
    for (unsigned int i = 0; i < current_iterate.multipliers.lower_bounds.size(); i++) {
        trial_multipliers.lower_bounds[i] = current_iterate.multipliers.lower_bounds[i] + dual_length * lower_delta_z[i];
        trial_multipliers.upper_bounds[i] = current_iterate.multipliers.upper_bounds[i] + dual_length * upper_delta_z[i];
        // TODO rescale the multipliers (IPOPT paper p6)
    }

    DEBUG << "MA57 restoration solution:\n";
    DEBUG << "Δx: "; print_vector(DEBUG, solution_IPM, 0, problem.number_variables);
    DEBUG << "Δs: "; print_vector(DEBUG, solution_IPM, problem.number_variables, problem.inequality_constraints.size());
    DEBUG << "Δz_L: "; print_vector(DEBUG, lower_delta_z);
    DEBUG << "Δz_U: "; print_vector(DEBUG, upper_delta_z);
    DEBUG << "primal length = " << primal_length << "\n";
    DEBUG << "dual length = " << dual_length << "\n\n";

    SubproblemSolution solution(trial_x, trial_multipliers);
    solution.status = INFEASIBLE;
    solution.phase = RESTORATION;
    solution.norm = norm_inf(solution.x, problem.number_variables);
    solution.predicted_reduction = [&](double step_length) {
        return this->compute_predicted_reduction(solution, step_length);
    };
    return solution;
}

bool InteriorPoint::phase_1_required(SubproblemSolution& solution) {
    return solution.status == INFEASIBLE;
}

double InteriorPoint::compute_central_complementarity_error(Iterate& iterate, double mu, std::vector<Range>& variables_bounds) {
    std::vector<double> residuals(iterate.x.size());
    /* variable bound constraints */
    for (unsigned int i = 0; i < iterate.x.size(); i++) {
        if (-INFINITY < variables_bounds[i].lb) {
            residuals[i] = iterate.multipliers.lower_bounds[i] * (iterate.x[i] - variables_bounds[i].lb) - mu;
        }
        if (variables_bounds[i].ub < INFINITY) {
            residuals[i] = iterate.multipliers.upper_bounds[i] * (iterate.x[i] - variables_bounds[i].ub) - mu;
        }
    }

    /* scaling */
    double sc = std::max(this->parameters.smax, (norm_1(iterate.multipliers.lower_bounds) + norm_1(iterate.multipliers.upper_bounds)) / iterate.x.size()) / this->parameters.smax;
    return norm(residuals, this->residual_norm) / sc;
}