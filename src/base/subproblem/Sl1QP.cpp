#include <cmath>
#include <map>
#include <memory>
#include "Sl1QP.hpp"
#include "Constraint.hpp"
#include "Utils.hpp"
#include "BQPDSolver.hpp"
#include "QPSolverFactory.hpp"

Sl1QP::Sl1QP(Problem& problem, std::string QP_solver, std::string hessian_evaluation_method, bool use_trust_region, bool scale_residuals):
// compute the number of variables and call the private constructor
Sl1QP(problem, QP_solver, hessian_evaluation_method, use_trust_region, scale_residuals, this->count_additional_variables_(problem)) {
}

int Sl1QP::count_additional_variables_(Problem& problem) {
    int number_variables = problem.number_variables;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (-INFINITY < problem.constraint_bounds[j].lb) {
            number_variables++;
        }
        if (problem.constraint_bounds[j].ub < INFINITY) {
            number_variables++;
        }
    }
    return number_variables;
}

Sl1QP::Sl1QP(Problem& problem, std::string QP_solver, std::string hessian_evaluation_method, bool use_trust_region, bool scale_residuals, int number_variables):
ActiveSetMethod(problem, QPSolverFactory::create(QP_solver, number_variables, problem.number_constraints, problem.hessian_maximum_number_nonzeros + problem.number_variables, true), scale_residuals),
// maximum number of Hessian nonzeros = number nonzeros + possible diagonal inertia correction
hessian_evaluation(HessianEvaluationFactory::create(hessian_evaluation_method, problem.number_variables)),
penalty_parameter(1.),
parameters({10., 0.1, 0.1}),
number_variables(number_variables) {
    // p and n are generated on the fly to solve the QP
    int current_index = problem.number_variables;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (-INFINITY < problem.constraint_bounds[j].lb) {
            // nonnegative variable p that captures the positive part of the constraint violation
            this->negative_part_variables[j] = current_index;
            current_index++;
        }
        if (problem.constraint_bounds[j].ub < INFINITY) {
            // nonnegative variable p that captures the negative part of the constraint violation
            this->positive_part_variables[j] = current_index;
            current_index++;
        }
    }
    
    /* if no trust region is used, the problem should be convexified by changing the inertia of the Hessian */
    this->hessian_evaluation->convexify = !use_trust_region;
}

SubproblemSolution Sl1QP::compute_step(Problem& problem, Iterate& current_iterate, double trust_region_radius) {
    DEBUG << "penalty parameter: " << this->penalty_parameter << "\n";

    // evaluate constraints
    current_iterate.compute_constraints(problem);

    /* stage a: compute the step within trust region */
    SubproblemSolution solution = this->solve_l1qp_subproblem_(problem, current_iterate, trust_region_radius, this->penalty_parameter);
    DEBUG << solution;

    /* penalty update: if penalty parameter is already 0, no need to decrease it */
    if (0. < this->penalty_parameter) {
        /* check infeasibility */
        double linearized_residual = this->compute_linearized_constraint_residual_(solution.x);
        DEBUG << "Linearized residual mk(dk): " << linearized_residual << "\n\n";

        // if problem had to be relaxed
        if (linearized_residual != 0.) {
            double current_penalty_parameter = this->penalty_parameter;

            /* stage c: solve the ideal l1 penalty problem with a zero penalty (no objective) */
            DEBUG << "Compute ideal solution:\n";
            SubproblemSolution ideal_solution = this->solve_l1qp_subproblem_(problem, current_iterate, trust_region_radius, 0.);
            DEBUG << ideal_solution;

            /* compute the ideal error (with a zero penalty parameter) */
            double ideal_error = this->compute_error_(problem, current_iterate, ideal_solution.multipliers, 0.);
            DEBUG << "Ideal error: " << ideal_error << "\n";

            if (ideal_error == 0.) {
                /* stage f: update the penalty parameter */
                this->penalty_parameter = 0.;
            }
            else {
                double ideal_linearized_residual = this->compute_linearized_constraint_residual_(ideal_solution.x);
                DEBUG << "Linearized residual mk(dk): " << ideal_linearized_residual << "\n\n";

                /* decrease penalty parameter to satisfy 2 conditions */
                bool condition1 = false, condition2 = false;
                while (!condition2) {
                    if (!condition1) {
                        /* stage d: reach a fraction of the ideal decrease */
                        if ((ideal_linearized_residual == 0. && linearized_residual == 0) || (ideal_linearized_residual != 0. && current_iterate.feasibility_measure - linearized_residual >= this->parameters.epsilon1 * (current_iterate.feasibility_measure - ideal_linearized_residual))) {
                            condition1 = true;
                            DEBUG << "Condition 1 is true\n";
                        }
                    }
                    /* stage e: further decrease penalty parameter if necessary */
                    if (condition1 && current_iterate.feasibility_measure - solution.objective >= this->parameters.epsilon2 * (current_iterate.feasibility_measure - ideal_solution.objective)) {
                        condition2 = true;
                        DEBUG << "Condition 2 is true\n";
                    }
                    if (!condition2) {
                        this->penalty_parameter /= this->parameters.tau;
                        if (this->penalty_parameter < 1e-10) {
                            this->penalty_parameter = 0.;
                            condition2 = true;
                        }
                        else {
                            DEBUG << "\nAttempting to solve with penalty parameter " << this->penalty_parameter << "\n";
                            solution = this->solve_l1qp_subproblem_(problem, current_iterate, trust_region_radius, this->penalty_parameter);
                            DEBUG << solution;

                            linearized_residual = this->compute_linearized_constraint_residual_(solution.x);
                            DEBUG << "Linearized residual mk(dk): " << linearized_residual << "\n\n";
                        }
                    }
                }

                /* stage f: update the penalty parameter */
                double term = ideal_error / std::max(1., current_iterate.feasibility_measure);
                this->penalty_parameter = std::min(this->penalty_parameter, term * term);
            }

            if (this->penalty_parameter < current_penalty_parameter) {
                DEBUG << "\n*** Penalty parameter updated to " << this->penalty_parameter << "\n";
                this->subproblem_definition_changed = true;
                /* recompute the solution */
                if (this->penalty_parameter == 0.) {
                    solution = ideal_solution;
                }
            }
        }
    }

    /* remove p and n */
    solution.x.resize(current_iterate.x.size());
    solution.multipliers.lower_bounds.resize(current_iterate.x.size());
    solution.multipliers.upper_bounds.resize(current_iterate.x.size());
    solution.norm = norm_inf(solution.x);

    /* remove contribution of positive part variables */
    for (std::pair<const int, int>& element: this->positive_part_variables) {
        int j = element.first;
        int i = element.second;
        current_iterate.constraints_jacobian[j].erase(i);
        current_iterate.objective_gradient.erase(i);
    }
    /* remove contribution of negative part variables */
    for (std::pair<const int, int>& element: this->negative_part_variables) {
        int j = element.first;
        int i = element.second;
        current_iterate.constraints_jacobian[j].erase(i);
        current_iterate.objective_gradient.erase(i);
    }
    return solution;
}

SubproblemSolution Sl1QP::solve_l1qp_subproblem_(Problem& problem, Iterate& current_iterate, double trust_region_radius, double penalty_parameter) {
    /* compute l1QP step */
    this->evaluate_optimality_iterate_(problem, current_iterate, penalty_parameter);
    SubproblemSolution solution = this->compute_qp_step_(problem, current_iterate, trust_region_radius);
    
    solution.objective_multiplier = penalty_parameter;
    solution.predicted_reduction = [&](double step_length) {
        return this->compute_predicted_reduction_(problem, current_iterate, solution, step_length);
    };
    
    // recompute active set: constraints are active when p-n = 0
    //this->recover_active_set(problem, solution, variables_bounds);
    return solution;
}

void Sl1QP::evaluate_optimality_iterate_(Problem& problem, Iterate& current_iterate, double penalty_parameter) {
    current_iterate.compute_objective_gradient(problem);
    std::map<int, double> objective_gradient;
    if (penalty_parameter != 0.) {
        for (std::pair<const int, double>& element: current_iterate.objective_gradient) {
            int i = element.first;
            double derivative = element.second;
            objective_gradient[i] = penalty_parameter*derivative;
        }
    }
    /* add contribution of positive part variables */
    for (std::pair<const int, int>& element: this->positive_part_variables) {
        int j = element.first;
        int i = element.second;
        current_iterate.constraints_jacobian[j][i] = -1.;
        objective_gradient[i] = 1.;
    }
    /* add contribution of negative part variables */
    for (std::pair<const int, int>& element: this->negative_part_variables) {
        int j = element.first;
        int i = element.second;
        current_iterate.constraints_jacobian[j][i] = 1.;
        objective_gradient[i] = 1.;
    }
    current_iterate.set_objective_gradient(objective_gradient);
    // Hessian
    current_iterate.is_hessian_computed = false;
    this->hessian_evaluation->compute(problem, current_iterate, penalty_parameter, current_iterate.multipliers.constraints);
    return;
}

SubproblemSolution Sl1QP::restore_feasibility(Problem&, Iterate&, SubproblemSolution&, double) {
    throw std::out_of_range("Sl1QP.compute_infeasibility_step is not implemented, since l1QP are always feasible");
}

double Sl1QP::compute_predicted_reduction_(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution, double step_length) {
    // the predicted reduction is quadratic
    if (step_length == 1.) {
        return current_iterate.feasibility_measure - solution.objective;
    }
    else {
        double linear_term = dot(solution.x, current_iterate.objective_gradient);
        double quadratic_term = current_iterate.hessian.quadratic_product(solution.x, solution.x) / 2.;
        // determine the constraint violation term: c(x_k) + alpha*\nabla c(x_k)^T d
        std::vector<double> scaled_constraints(current_iterate.constraints);
        for (unsigned int j = 0; j < current_iterate.constraints.size(); j++) {
            scaled_constraints[j] += step_length * dot(solution.x, current_iterate.constraints_jacobian[j]);
        }
        double constraint_violation = problem.compute_constraint_residual(scaled_constraints, this->residual_norm);
        return current_iterate.feasibility_measure - (step_length * (linear_term + step_length * quadratic_term) + constraint_violation);
    }
}

/* private methods */

std::vector<Range> Sl1QP::generate_variables_bounds_(Problem& problem, Iterate& current_iterate, double trust_region_radius) {
    // p and n are nonnegative
    std::vector<Range> variables_bounds(this->number_variables, {0., INFINITY});

    /* original bounds intersected with trust region  */
    for (int i = 0; i < problem.number_variables; i++) {
        double lb = std::max(-trust_region_radius, this->subproblem_variables_bounds[i].lb - current_iterate.x[i]);
        double ub = std::min(trust_region_radius, this->subproblem_variables_bounds[i].ub - current_iterate.x[i]);
        variables_bounds[i] = {lb, ub};
    }
    return variables_bounds;
}

double Sl1QP::compute_linearized_constraint_residual_(std::vector<double>& d) {
    double residual = 0.;
    // l1 residual of the linearized constraints
    for (std::pair<const int, int>& element: this->positive_part_variables) {
        int i = element.second;
        residual += d[i];
    }
    for (std::pair<const int, int>& element: this->negative_part_variables) {
        int i = element.second;
        residual += d[i];
    }
    return residual;
}

double Sl1QP::compute_error_(Problem& problem, Iterate& iterate, Multipliers& multipliers, double penalty_parameter) {
    /* measure that combines KKT error and complementarity error */
    double error = 0.;

    /* KKT error */
    std::vector<double> lagrangian_gradient = iterate.lagrangian_gradient(problem, penalty_parameter, multipliers);
    error += norm_1(lagrangian_gradient);
    /* complementarity error */
    error += this->compute_complementarity_error_(problem, iterate, multipliers);
    return error;
}

/* complementary slackness error. Use abs/1e-8 to safeguard */
double Sl1QP::compute_complementarity_error_(Problem& problem, Iterate& iterate, Multipliers& multipliers) {
    double error = 0.;
    /* bound constraints */
    for (int i = 0; i < problem.number_variables; i++) {
        if (-INFINITY < problem.variables_bounds[i].lb) {
            error += std::abs(iterate.multipliers.lower_bounds[i] * (iterate.x[i] - problem.variables_bounds[i].lb));
        }
        if (problem.variables_bounds[i].ub < INFINITY) {
            error += std::abs(iterate.multipliers.upper_bounds[i] * (iterate.x[i] - problem.variables_bounds[i].ub));
        }
    }
    /* general constraints */
    for (int j = 0; j < problem.number_constraints; j++) {
        double multiplier_j = multipliers.constraints[j];
        if (iterate.constraints[j] < problem.constraint_bounds[j].lb) {
            // violated lower: the multiplier is 1 at optimum
            error += std::abs((1. - multiplier_j) * (problem.constraint_bounds[j].lb - iterate.constraints[j]));
        }
        else if (problem.constraint_bounds[j].ub < iterate.constraints[j]) {
            // violated upper: the multiplier is -1 at optimum
            error += std::abs((1. + multiplier_j) * (iterate.constraints[j] - problem.constraint_bounds[j].ub));
        }
        else if (-INFINITY < problem.constraint_bounds[j].lb && 0. < multiplier_j) {
            error += std::abs(multiplier_j * (iterate.constraints[j] - problem.constraint_bounds[j].lb));
        }
        else if (problem.constraint_bounds[j].ub < INFINITY && multiplier_j < 0.) {
            error += std::abs(multiplier_j * (iterate.constraints[j] - problem.constraint_bounds[j].ub));
        }
    }
    return error;
}

void Sl1QP::recover_active_set_(Problem& problem, SubproblemSolution& solution, std::vector<Range>& variables_bounds) {
    solution.active_set.at_lower_bound.clear();
    solution.active_set.at_upper_bound.clear();
    // variables
    for (int i = 0; i < problem.number_variables; i++) {
        if (solution.x[i] == variables_bounds[i].lb) {
            solution.active_set.at_lower_bound.insert(i);
        }
        else if (solution.x[i] == variables_bounds[i].ub) {
            solution.active_set.at_upper_bound.insert(i);
        }
    }
    // constraints: only when p-n = 0
    for (unsigned int j = 0; j < solution.multipliers.constraints.size(); j++) {
        double constraint_violation = 0.;
        if (positive_part_variables.find(j) != positive_part_variables.end()) {
            constraint_violation += solution.x[this->positive_part_variables[j]];
        }
        if (negative_part_variables.find(j) != negative_part_variables.end()) {
            constraint_violation += solution.x[this->negative_part_variables[j]];
        }
        if (constraint_violation == 0.) {
            solution.active_set.at_lower_bound.insert(problem.number_variables + j);
        }
    }
    return;
}