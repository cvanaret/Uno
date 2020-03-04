#include <cmath>
#include "Argonot.hpp"
#include "PenaltyMeritFunction.hpp"
#include "Logger.hpp"

/*
 * Infeasibility detection and SQP methods for nonlinear optimization 
 * http://epubs.siam.org/doi/pdf/10.1137/080738222
 */

PenaltyMeritFunction::PenaltyMeritFunction(Subproblem& subproblem, double tolerance) :
GlobalizationStrategy(subproblem, tolerance), penalty_parameter(1.) {
    this->tau = 0.5;
    this->eta = 1e-8;
    this->epsilon1 = 0.1;
    this->epsilon2 = 0.1;
}

Iterate PenaltyMeritFunction::initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers, bool use_trust_region) {
    /* compute the number of necessary additional variables and constraints */
    this->penalty_dimensions.number_additional_variables = 0;
    this->penalty_dimensions.number_constraints = 0;

    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] == EQUAL_BOUNDS) {
            this->penalty_dimensions.number_additional_variables += 2;
            this->penalty_dimensions.number_constraints++;
        } else if (problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
            this->penalty_dimensions.number_additional_variables += 2;
            this->penalty_dimensions.number_constraints += 2;
        } else if (problem.constraint_status[j] == BOUNDED_LOWER) {
            this->penalty_dimensions.number_additional_variables++;
            this->penalty_dimensions.number_constraints++;
        } else {
            this->penalty_dimensions.number_additional_variables++;
            this->penalty_dimensions.number_constraints++;
        }
    }

    /* allocate the subproblem solver */
    int number_variables = problem.number_variables + this->penalty_dimensions.number_additional_variables;
    int number_constraints = this->penalty_dimensions.number_constraints;
    Iterate first_iterate = this->subproblem.initialize(problem, x, multipliers, number_variables, number_constraints, use_trust_region);

    first_iterate.KKTerror = Argonot::compute_KKT_error(problem, first_iterate, this->penalty_parameter);
    first_iterate.complementarity_error = Argonot::compute_complementarity_error(problem, first_iterate);

    return first_iterate;
}

SubproblemSolution PenaltyMeritFunction::compute_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds) {
    /* stage a: compute the step within trust region */
    SubproblemSolution solution = this->subproblem.compute_l1_penalty_step(problem, current_iterate, variables_bounds, this->penalty_parameter, this->penalty_dimensions);
    DEBUG << solution;

    /* if penalty parameter is already 0, no need to decrease it */
    if (0. < this->penalty_parameter) {
        /* check infeasibility */
        double linear_model = this->compute_linear_model(problem, solution);
        if (linear_model != 0.) {
            double current_penalty_parameter = this->penalty_parameter;

            /* stage c: solve the ideal l1 penalty problem with a zero penalty (no objective) */
            SubproblemSolution ideal_solution = this->subproblem.compute_l1_penalty_step(problem, current_iterate, variables_bounds, 0., this->penalty_dimensions);
            DEBUG << ideal_solution;

            /* stage f: update the penalty parameter */
            std::vector<double> ideal_bound_multipliers = this->compute_bound_multipliers(problem, ideal_solution);
            std::vector<double> ideal_constraint_multipliers = this->compute_constraint_multipliers(problem, ideal_solution);
            Multipliers ideal_multipliers = {ideal_bound_multipliers, ideal_constraint_multipliers};
            /* compute the ideal error (with a zero penalty parameter) */
            double ideal_error = this->compute_error(problem, current_iterate, ideal_multipliers, 0.);

            if (ideal_error == 0.) {
                /* stage f: update the penalty parameter */
                this->penalty_parameter = 0.;
            } else {
                double ideal_linear_model = this->compute_linear_model(problem, ideal_solution);

                /* decrease penalty parameter to satisfy 2 conditions */
                bool condition1 = false, condition2 = false;
                while (!condition2) {
                    this->penalty_parameter *= this->tau;
                    if (this->penalty_parameter < 1e-10) {
                        this->penalty_parameter = 0.;
                        condition2 = true;
                    }

                    DEBUG << "\nSolving with penalty parameter " << this->penalty_parameter << "\n";
                    solution = this->subproblem.compute_l1_penalty_step(problem, current_iterate, variables_bounds, this->penalty_parameter, this->penalty_dimensions);
                    DEBUG << solution;

                    double trial_linear_model = this->compute_linear_model(problem, solution);
                    if (!condition1) {
                        /* stage d: reach a fraction of the ideal decrease */
                        if ((ideal_linear_model == 0. && trial_linear_model == 0.) || (ideal_linear_model != 0. &&
                                current_iterate.residual - trial_linear_model >= this->epsilon1 * (current_iterate.residual - ideal_linear_model))) {
                            condition1 = true;
                        }
                    }
                    /* stage e: further decrease penalty parameter if necessary */
                    if (condition1 && current_iterate.residual - solution.objective >= this->epsilon2 * (current_iterate.residual - ideal_solution.objective)) {
                        condition2 = true;
                    }
                }

                /* stage f: update the penalty parameter */
                double term = ideal_error / std::max(1., current_iterate.residual);
                this->penalty_parameter = std::min(this->penalty_parameter, term * term);
            }

            if (this->penalty_parameter < current_penalty_parameter) {
                DEBUG << "Penalty parameter updated to " << this->penalty_parameter << "\n";
                /* recompute the solution */
                if (this->penalty_parameter == 0.) {
                    solution = ideal_solution;
                } else {
                    solution = this->subproblem.compute_l1_penalty_step(problem, current_iterate, variables_bounds, this->penalty_parameter, this->penalty_dimensions);
                    DEBUG << solution;
                }
            }
        }
    }
    INFO << "penalty parameter: " << this->penalty_parameter << "\t";
    return solution;
}

bool PenaltyMeritFunction::check_step(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution, double step_length) {
    /* stage g: line-search along fixed step */

    /* retrieve only original primal and dual variables from the step */
    std::vector<double> d(problem.number_variables);
    for (int i = 0; i < problem.number_variables; i++) {
        d[i] = solution.x[i];
    }

    /* generate the trial point */
    std::vector<double> x_trial = add_vectors(current_iterate.x, d, step_length);
    /* get the multipliers */
    std::vector<double> bound_multipliers = this->compute_bound_multipliers(problem, solution);
    std::vector<double> constraint_multipliers = this->compute_constraint_multipliers(problem, solution);
    std::vector<double> trial_constraint_multipliers = add_vectors(current_iterate.multipliers.constraints, constraint_multipliers, step_length);
    Multipliers trial_multipliers = {bound_multipliers, trial_constraint_multipliers};
    /* generate the trial iterate */
    Iterate trial_iterate(problem, x_trial, trial_multipliers);
    this->compute_measures(problem, trial_iterate);

    /* compute current exact l1 penalty: rho f + sum max(0, c) */
    double current_exact_l1_penalty = this->penalty_parameter * current_iterate.objective + current_iterate.residual;
    /* compute trial exact l1 penalty */
    double trial_exact_l1_penalty = this->penalty_parameter * trial_iterate.objective + trial_iterate.residual;

    /* check the validity of the trial step */
    bool accept = false;
    if (current_exact_l1_penalty - trial_exact_l1_penalty >= this->eta * step_length * (current_iterate.residual - solution.objective)) {
        accept = true;
        trial_iterate.KKTerror = Argonot::compute_KKT_error(problem, trial_iterate, this->penalty_parameter);
        trial_iterate.complementarity_error = Argonot::compute_complementarity_error(problem, trial_iterate);
        double step_norm = step_length * norm_inf(d);
        trial_iterate.status = this->compute_status(problem, trial_iterate, step_norm);
        current_iterate = trial_iterate;
    }
    return accept;
}

OptimalityStatus PenaltyMeritFunction::compute_status(Problem& problem, Iterate& trial_iterate, double step_norm) {
    OptimalityStatus status = NOT_OPTIMAL;

    /* test for optimality */
    double optimality_error = this->compute_error(problem, trial_iterate, trial_iterate.multipliers, this->penalty_parameter);
    DEBUG << "Ek(lambda_k, rho_k) = " << optimality_error << "\n";
    if (optimality_error <= this->tolerance && trial_iterate.residual <= this->tolerance * problem.number_constraints) {
        status = KKT_POINT;
        /* rescale the multipliers */
        if (0. < this->penalty_parameter) {
            for (int i = 0; i < problem.number_variables; i++) {
                trial_iterate.multipliers.lower_bounds[i] /= this->penalty_parameter;
                trial_iterate.multipliers.upper_bounds[i] /= this->penalty_parameter;
            }
            for (int j = 0; j < problem.number_constraints; j++) {
                trial_iterate.multipliers.constraints[j] /= this->penalty_parameter;
            }
        }
    }
    else {
        double infeasibility_error = this->compute_error(problem, trial_iterate, trial_iterate.multipliers, 0.);
        DEBUG << "Ek(lambda_k, 0.) = " << infeasibility_error << "\n";
        if (infeasibility_error <= this->tolerance && trial_iterate.residual > this->tolerance * problem.number_constraints) {
            status = FJ_POINT;
        } else if (step_norm <= this->tolerance / 100.) {
            if (trial_iterate.residual <= this->tolerance * problem.number_constraints) {
                status = FEASIBLE_SMALL_STEP;
            } else {
                status = INFEASIBLE_SMALL_STEP;
            }
        }
    }
    return status;
}

double PenaltyMeritFunction::compute_linear_model(Problem& problem, SubproblemSolution& solution) {
    double linear_model = 0.;
    for (int k = 0; k < this->penalty_dimensions.number_additional_variables; k++) {
        linear_model += solution.x[problem.number_variables + k];
    }
    return linear_model;
}

std::vector<double> PenaltyMeritFunction::compute_bound_multipliers(Problem& problem, SubproblemSolution& solution) {
    std::vector<double> bound_multipliers(problem.number_variables);
    for (int i = 0; i < problem.number_variables; i++) {
        bound_multipliers[i] = solution.multipliers.lower_bounds[i] + solution.multipliers.upper_bounds[i];
    }
    return bound_multipliers;
}

std::vector<double> PenaltyMeritFunction::compute_constraint_multipliers(Problem& problem, SubproblemSolution& solution) {
    std::vector<double> constraint_multipliers(problem.number_constraints);
    int current_constraint = 0;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] == BOUNDED_BOTH_SIDES) {
            /* only case where 2 constraints were generated */
            /* only one bound is active: one multiplier is > 0, the other is 0 */
            constraint_multipliers[j] = solution.multipliers.constraints[current_constraint] + solution.multipliers.constraints[current_constraint + 1];
            current_constraint += 2;
        } else {
            /* only 1 constraint was generated */
            constraint_multipliers[j] = solution.multipliers.constraints[current_constraint];
            current_constraint++;
        }
    }
    return constraint_multipliers;
}

double PenaltyMeritFunction::compute_error(Problem& problem, Iterate& current_iterate, Multipliers& multipliers, double penalty_parameter) {
    /* measure that combines KKT error and complementarity error */
    double error = 0.;

    /* KKT error */
    std::vector<double> lagrangian_gradient = Argonot::compute_lagrangian_gradient(problem, current_iterate, penalty_parameter, multipliers);
    // compute 1-norm
    error += norm_1(lagrangian_gradient);

    /* complementarity error */
    // bound constraints
    for (int i = 0; i < problem.number_variables; i++) {
        if (problem.variables_bounds[i].lb < current_iterate.x[i] && current_iterate.x[i] < problem.variables_bounds[i].ub) {
            error += std::abs(multipliers.lower_bounds[i] * (current_iterate.x[i] - problem.variables_bounds[i].lb));
            error += std::abs(multipliers.upper_bounds[i] * (current_iterate.x[i] - problem.variables_bounds[i].ub));
        }
    }
    // check if constraint is strictly satisfied or violated
    for (int j = 0; j < problem.number_constraints; j++) {
        double multiplier_j = multipliers.constraints[j];

        /* violated */
        if (current_iterate.constraints[j] < problem.constraints_bounds[j].lb) {
            error += std::abs((1. - multiplier_j)*(current_iterate.constraints[j] - problem.constraints_bounds[j].lb));
        } else if (problem.constraints_bounds[j].ub < current_iterate.constraints[j]) {
            error += std::abs((1. - multiplier_j)*(current_iterate.constraints[j] - problem.constraints_bounds[j].ub));

        } else {
            /* active or strictly satisfied */
            if (multiplier_j > 0.) {
                error += std::abs(multiplier_j * (current_iterate.constraints[j] - problem.constraints_bounds[j].lb));
            } else if (multiplier_j < 0.) {
                error += std::abs(multiplier_j * (current_iterate.constraints[j] - problem.constraints_bounds[j].ub));
            }
        }
    }
    return error;
}

//for (int j = 0; j < problem.number_constraints; j++) {
//double multiplier_j = multipliers[problem.number_variables + j];
///* strictly satisfied */
//if (problem.constraints_bounds[j].lb < current_iterate.constraints[j] && current_iterate.constraints[j] < problem.constraints_bounds[j].ub) {
//if (multiplier_j > 0.) {
//error += std::abs(multiplier_j*(current_iterate.constraints[j] - problem.constraints_bounds[j].lb));
//}
//else if (multiplier_j < 0.) {
//error += std::abs(multiplier_j*(current_iterate.constraints[j] - problem.constraints_bounds[j].ub));
//}
//}
///* violated */
//else if (current_iterate.constraints[j] < problem.constraints_bounds[j].lb) {
//error += std::abs((1. - multiplier_j)*(current_iterate.constraints[j] - problem.constraints_bounds[j].lb));
//}
//else if (problem.constraints_bounds[j].ub < current_iterate.constraints[j]) {
//error += std::abs((1. - multiplier_j)*(current_iterate.constraints[j] - problem.constraints_bounds[j].ub));

//}
//}

void PenaltyMeritFunction::compute_measures(Problem& problem, Iterate& iterate) {
    this->subproblem.compute_measures(problem, iterate);
    // TODO: add the penalized term to the optimality measure
    return;
}
