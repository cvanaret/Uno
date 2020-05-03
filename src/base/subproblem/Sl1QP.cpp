#include <cmath>
#include <map>
#include "Sl1QP.hpp"
#include "Constraint.hpp"
#include "Utils.hpp"

Sl1QP::Sl1QP(QPSolver& solver, HessianEvaluation& hessian_evaluation): Subproblem("l1"), solver(solver), hessian_evaluation(hessian_evaluation), penalty_parameter(1.), parameters({10., 0.1, 0.1}) {
}

Iterate Sl1QP::initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers, bool use_trust_region) {
    // register the original bounds
    this->subproblem_variables_bounds = problem.variables_bounds;
    
    // p and n are generated on the fly to solve the QP, but are not kept
    int current_index = problem.number_variables;
    for (int j = 0; j < problem.number_constraints; j++) {
        // nonnegative variable p that captures the positive part of an equality
        this->positive_part_variables[j] = current_index;
        current_index++;
        // nonnegative variable p that captures the positive part of an equality
        this->negative_part_variables[j] = current_index;
        current_index++;
    }

    /* compute the optimality and feasibility measures of the initial point */
    Iterate first_iterate(x, multipliers);
    this->compute_optimality_measures(problem, first_iterate);

    /* allocate the QP solver */
    // p and n are generated on the fly to solve the QP, but are not kept
    //int number_variables_qp = problem.number_variables + 2 * problem.number_constraints;
    //this->solver.allocate(number_variables_qp, problem.number_constraints);
    
    /* if no trust region is used, the problem should be convexified by changing the inertia of the Hessian */
    this->hessian_evaluation.convexify = !use_trust_region;
    
    return first_iterate;
}

SubproblemSolution Sl1QP::compute_optimality_step(Problem& problem, Iterate& current_iterate, double trust_region_radius) {
    DEBUG << "penalty parameter: " << this->penalty_parameter << "\n";

    // evaluate constraints
    current_iterate.compute_constraints(problem);

    /* stage a: compute the step within trust region */
    SubproblemSolution solution = this->solve_subproblem(problem, current_iterate, trust_region_radius, this->penalty_parameter);
    DEBUG << solution;

    /* penalty update: if penalty parameter is already 0, no need to decrease it */
    if (0. < this->penalty_parameter) {
        /* check infeasibility */
        double linearized_residual = this->compute_linearized_constraint_residual(problem, solution.x);
        DEBUG << "Linearized residual mk(dk): " << linearized_residual << "\n\n";

        // if problem had to be relaxed
        if (linearized_residual != 0.) {
            double current_penalty_parameter = this->penalty_parameter;

            /* stage c: solve the ideal l1 penalty problem with a zero penalty (no objective) */
            DEBUG << "Compute ideal solution:\n";
            SubproblemSolution ideal_solution = this->solve_subproblem(problem, current_iterate, trust_region_radius, 0.);
            DEBUG << ideal_solution;

            /* compute the ideal error (with a zero penalty parameter) */
            double ideal_error = this->compute_error(problem, current_iterate, ideal_solution.multipliers, 0.);
            DEBUG << "Ideal error: " << ideal_error << "\n";

//            /* stage f: update the penalty parameter */
//            double term = ideal_error / std::max(1., current_iterate.feasibility_measure);
//            this->penalty_parameter = std::min(this->penalty_parameter, term * term);
            
            if (ideal_error == 0.) {
                /* stage f: update the penalty parameter */
                this->penalty_parameter = 0.;
            }
            else {
                double ideal_linearized_residual = this->compute_linearized_constraint_residual(problem, ideal_solution.x);
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
                            solution = this->solve_subproblem(problem, current_iterate, trust_region_radius, this->penalty_parameter);
                            DEBUG << solution;

                            linearized_residual = this->compute_linearized_constraint_residual(problem, solution.x);
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
    //INFO << "penalty parameter: " << this->penalty_parameter << "\t";

    /* remove p and n */
    std::vector<double> filtered_x(solution.x.begin(), solution.x.begin() + current_iterate.x.size());
    std::vector<double> filtered_lb_multipliers(solution.multipliers.lower_bounds.begin(), solution.multipliers.lower_bounds.begin() + current_iterate.x.size());
    std::vector<double> filtered_ub_multipliers(solution.multipliers.upper_bounds.begin(), solution.multipliers.upper_bounds.begin() + current_iterate.x.size());
    solution.x = filtered_x;
    solution.norm = norm_inf(filtered_x);
    solution.multipliers.lower_bounds = filtered_lb_multipliers;
    solution.multipliers.upper_bounds = filtered_ub_multipliers;
    
    /* remove contribution of positive part variables */
    for (std::pair<const int, int>& element: this->positive_part_variables) {
        int j = element.first;
        int i = element.second;
        current_iterate.constraints_jacobian[j].erase(i);
    }
    /* remove contribution of negative part variables */
    for (std::pair<const int, int>& element: this->negative_part_variables) {
        int j = element.first;
        int i = element.second;
        current_iterate.constraints_jacobian[j].erase(i);
    }
    return solution;
}

SubproblemSolution Sl1QP::solve_subproblem(Problem& problem, Iterate& current_iterate, double trust_region_radius, double penalty_parameter) {
    DEBUG << "Current point: "; print_vector(DEBUG, current_iterate.x);

    /* compute first- and second-order information */
    current_iterate.compute_objective_gradient(problem);
    std::map<int, double> objective_gradient;
    if (penalty_parameter != 0.) {
        for (std::pair<const int, double>& element: current_iterate.objective_gradient) {
            int i = element.first;
            int derivative = element.second;
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
    // Hessian
    current_iterate.is_hessian_computed = false;
    this->hessian_evaluation.compute(problem, current_iterate, penalty_parameter, current_iterate.multipliers.constraints);

    /* bounds of the variables */
    std::vector<Range> variables_bounds = this->generate_variables_bounds(problem, current_iterate, trust_region_radius);

    /* bounds of the linearized equality constraints */
    std::vector<Range> constraints_bounds = Subproblem::generate_constraints_bounds(problem, current_iterate.constraints);
    
    DEBUG << "hessian: " << current_iterate.hessian;
    DEBUG << "gradient obj: ";
    print_vector(DEBUG, objective_gradient);
    for (int j = 0; j < problem.number_constraints; j++) {
        DEBUG << "gradient c" << j << ": ";
        print_vector(DEBUG, current_iterate.constraints_jacobian[j]);
    }
    for (unsigned int i = 0; i < current_iterate.x.size(); i++) {
        DEBUG << "x" << i << " in [" << this->subproblem_variables_bounds[i].lb << ", " << this->subproblem_variables_bounds[i].ub << "]\n";
    }
    for (unsigned int i = 0; i < variables_bounds.size(); i++) {
        DEBUG << "Δx" << i << " in [" << variables_bounds[i].lb << ", " << variables_bounds[i].ub << "]\n";
    }
    for (unsigned int j = 0; j < constraints_bounds.size(); j++) {
        DEBUG << "linearized c" << j << " in [" << constraints_bounds[j].lb << ", " << constraints_bounds[j].ub << "]\n";
    }

    /* generate the initial point */
    std::vector<double> d0(variables_bounds.size()); // = {0.}

    /* solve the QP */
    SubproblemSolution solution = this->solver.solve_QP(variables_bounds, constraints_bounds, objective_gradient, current_iterate.constraints_jacobian, current_iterate.hessian, d0);

    // recompute active set: constraints are active when p-n = 0
    this->recover_active_set(problem, solution, variables_bounds);

    solution.phase_1_required = this->phase_1_required(solution);
    solution.objective_multiplier = penalty_parameter;
    this->number_subproblems_solved++;
    
    return solution;
}

SubproblemSolution Sl1QP::compute_infeasibility_step(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution, double trust_region_radius) {
    throw std::out_of_range("Sl1QP.compute_infeasibility_step is not implemented, since l1QP are always feasible");
}

double Sl1QP::compute_predicted_reduction(Iterate& current_iterate, SubproblemSolution& solution) {
    return current_iterate.feasibility_measure - solution.objective;
}

void Sl1QP::compute_optimality_measures(Problem& problem, Iterate& iterate) {
    /* feasibility */
    iterate.compute_constraint_residual(problem, this->residual_norm);
    iterate.feasibility_measure = iterate.constraint_residual;
    /* optimality */
    iterate.compute_objective(problem);
    iterate.optimality_measure = iterate.objective;
    return;
}

void Sl1QP::compute_infeasibility_measures(Problem& problem, Iterate& iterate, SubproblemSolution& solution) {
    throw std::out_of_range("Sl1QP.compute_infeasibility_measures is not implemented, since l1QP are always feasible");
}

bool Sl1QP::phase_1_required(SubproblemSolution& solution) {
    //TODO
    return false;
}

/* private methods */

std::vector<Range> Sl1QP::generate_variables_bounds(Problem& problem, Iterate& current_iterate, double trust_region_radius) {
    std::vector<Range> variables_bounds(current_iterate.x.size() + 2 * problem.number_constraints);
    /* original bounds intersected with trust region  */
    for (unsigned int i = 0; i < current_iterate.x.size(); i++) {
        double lb = std::max(-trust_region_radius, this->subproblem_variables_bounds[i].lb - current_iterate.x[i]);
        double ub = std::min(trust_region_radius, this->subproblem_variables_bounds[i].ub - current_iterate.x[i]);
        variables_bounds[i] = {lb, ub};
    }
    /* p and n are non-negative */
    for (unsigned int i = current_iterate.x.size(); i < current_iterate.x.size() + 2 * problem.number_constraints; i++) {
        variables_bounds[i] = {0., INFINITY};
    }
    return variables_bounds;
}

double Sl1QP::compute_linearized_constraint_residual(Problem& problem, std::vector<double>& x) {
    double residual = 0.;
    // l1 residual of the linearized constraints
    for (int j = 0; j < problem.number_constraints; j++) {
        residual += x[this->positive_part_variables[j]] + x[this->negative_part_variables[j]];
    }
    return residual;
}

double Sl1QP::compute_error(Problem& problem, Iterate& current_iterate, Multipliers& multipliers, double penalty_parameter) {
    /* measure that combines KKT error and complementarity error */
    double error = 0.;

    /* KKT error */
    std::vector<double> lagrangian_gradient = current_iterate.lagrangian_gradient(problem, penalty_parameter, multipliers);
    // compute 1-norm
    error += norm_1(lagrangian_gradient);

    /* complementarity error of bound constraints */
    for (int j = 0; j < problem.number_constraints; j++) {
        if (current_iterate.constraints[j] < problem.constraints_bounds[j].lb) {
            // violated lower: the multiplier is 1 at optimum
            error += std::abs((1. - multipliers.constraints[j]) * (problem.constraints_bounds[j].lb - current_iterate.constraints[j]));
        }
        else if (problem.constraints_bounds[j].ub < current_iterate.constraints[j]) {
            // violated upper: the multiplier is -1 at optimum
            error += std::abs((1. + multipliers.constraints[j]) * (current_iterate.constraints[j] - problem.constraints_bounds[j].ub));
        }
        else {
            // strictly satisfied: the multiplier is 0 at optimum
            error += std::abs(multipliers.constraints[j]*current_iterate.constraints[j]);
        }
    }
    return error;
}

void Sl1QP::recover_active_set(Problem& problem, SubproblemSolution& solution, std::vector<Range>& variables_bounds) {
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
        if (solution.x[this->positive_part_variables[j]] + solution.x[this->negative_part_variables[j]] == 0.) {
            solution.active_set.at_lower_bound.insert(problem.number_variables + j);
        }
    }
    return;
}