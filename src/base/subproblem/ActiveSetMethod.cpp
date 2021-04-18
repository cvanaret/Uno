#include <cmath>
#include <map>
#include "ActiveSetMethod.hpp"
#include "Constraint.hpp"
#include "Utils.hpp"
#include "Logger.hpp"

ActiveSetMethod::ActiveSetMethod(Problem& problem, bool scale_residuals):
Subproblem("l1", problem.variables_bounds, scale_residuals) {
}

Iterate ActiveSetMethod::evaluate_initial_point(Problem& problem, std::vector<double>& x, Multipliers& multipliers) {
    Iterate first_iterate(x, multipliers);
    /* compute the optimality and feasibility measures of the initial point */
    this->compute_optimality_measures(problem, first_iterate);
    return first_iterate;
}

std::vector<Range> ActiveSetMethod::generate_variables_bounds_(Problem& /*problem*/, Iterate& current_iterate, double trust_region_radius) {
    std::vector<Range> bounds(current_iterate.x.size());
    /* bounds intersected with trust region  */
    for (unsigned int i = 0; i < current_iterate.x.size(); i++) {
        double lb = std::max(-trust_region_radius, this->subproblem_variables_bounds[i].lb - current_iterate.x[i]);
        double ub = std::min(trust_region_radius, this->subproblem_variables_bounds[i].ub - current_iterate.x[i]);
        bounds[i] = {lb, ub};
    }
    return bounds;
}

void ActiveSetMethod::compute_optimality_measures(Problem& problem, Iterate& iterate) {
    // feasibility
    this->compute_residuals(problem, iterate, iterate.multipliers, 1.);
    iterate.feasibility_measure = iterate.residuals.constraints;
    // optimality
    iterate.compute_objective(problem);
    iterate.optimality_measure = iterate.objective;
    return;
}

void ActiveSetMethod::compute_infeasibility_measures(Problem& problem, Iterate& iterate, Direction& direction) {
    iterate.compute_constraints(problem);
    // feasibility measure: residual of all constraints
    iterate.feasibility_measure = problem.compute_constraint_residual(iterate.constraints, this->residual_norm);
    // optimality measure: residual of linearly infeasible constraints
    iterate.optimality_measure = problem.compute_constraint_residual(iterate.constraints, direction.constraint_partition.infeasible, this->residual_norm);
    return;
}

/* QP */

Direction ActiveSetMethod::compute_qp_step_(Problem& problem, QPSolver& solver, Iterate& current_iterate, double trust_region_radius) {
    DEBUG << "Current point: "; print_vector(DEBUG, current_iterate.x);
    DEBUG << "Current constraint multipliers: "; print_vector(DEBUG, current_iterate.multipliers.constraints);
    DEBUG << "Current lb multipliers: "; print_vector(DEBUG, current_iterate.multipliers.lower_bounds);
    DEBUG << "Current ub multipliers: "; print_vector(DEBUG, current_iterate.multipliers.upper_bounds);
    
    /* bounds of the variables */
    std::vector<Range> variables_bounds = this->generate_variables_bounds_(problem, current_iterate, trust_region_radius);

    /* bounds of the linearized constraints */
    std::vector<Range> constraints_bounds = Subproblem::generate_constraints_bounds(problem, current_iterate.constraints);

    /* generate the initial point */
    std::vector<double> d0(variables_bounds.size()); // = {0.}

    /* solve the QP */
    Direction direction = solver.solve_QP(variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate.constraints_jacobian, current_iterate.hessian, d0);
    this->number_subproblems_solved++;
    DEBUG << direction;
    return direction;
}

double ActiveSetMethod::compute_qp_predicted_reduction_(Problem& /*problem*/, Iterate& current_iterate, Direction& direction, double step_length) {
    // the predicted reduction is quadratic in the step length
    if (step_length == 1.) {
        return -direction.objective;
    }
    else {
        double linear_term = dot(direction.x, current_iterate.objective_gradient);
        double quadratic_term = current_iterate.hessian.quadratic_product(direction.x, direction.x) / 2.;
        return -step_length * (linear_term + step_length * quadratic_term);
    }
}

Direction ActiveSetMethod::compute_l1qp_step_(Problem& problem, QPSolver& solver, Iterate& current_iterate, ConstraintPartition& constraint_partition, std::vector<double>& initial_solution, double trust_region_radius) {
    /* compute the objective */
    this->compute_l1_linear_objective_(current_iterate, constraint_partition);

    /* bounds of the variables */
    std::vector<Range> variables_bounds = this->generate_variables_bounds_(problem, current_iterate, trust_region_radius);

    /* bounds of the linearized constraints */
    std::vector<Range> constraints_bounds = this->generate_feasibility_bounds_(problem, current_iterate.constraints, constraint_partition);

    /* generate the initial point */
    std::vector<double>& d0 = initial_solution;

    /* solve the QP */
    Direction direction = solver.solve_QP(variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate.constraints_jacobian, current_iterate.hessian, d0);
    direction.objective_multiplier = 0.;
    direction.constraint_partition = constraint_partition;
    this->number_subproblems_solved++;
    DEBUG << direction;
    return direction;
}

Direction ActiveSetMethod::compute_l1qp_step_(Problem& problem, QPSolver& solver, Iterate& current_iterate, double penalty_parameter, ElasticVariables& elastic_variables, double trust_region_radius) {
    current_iterate.compute_objective_gradient(problem);
    SparseGradient objective_gradient;
    if (penalty_parameter != 0.) {
        for (std::pair<const int, double>& element: current_iterate.objective_gradient) {
            int i = element.first;
            double derivative = element.second;
            objective_gradient[i] = penalty_parameter*derivative;
        }
    }
    /* add contribution of positive part variables */
    for (std::pair<const int, int>& element: elastic_variables.positive) {
        int j = element.first;
        int i = element.second;
        current_iterate.constraints_jacobian[j][i] = -1.;
        objective_gradient[i] = 1.;
    }
    /* add contribution of negative part variables */
    for (std::pair<const int, int>& element: elastic_variables.negative) {
        int j = element.first;
        int i = element.second;
        current_iterate.constraints_jacobian[j][i] = 1.;
        objective_gradient[i] = 1.;
    }
    current_iterate.set_objective_gradient(objective_gradient);
    
    /* bounds of the variables */
    std::vector<Range> variables_bounds = this->generate_variables_bounds_(problem, current_iterate, trust_region_radius);

    /* bounds of the linearized constraints */
    std::vector<Range> constraints_bounds = Subproblem::generate_constraints_bounds(problem, current_iterate.constraints);

    /* generate the initial point */
    std::vector<double> d0(variables_bounds.size()); // = {0.}

    /* solve the QP */
    Direction direction = solver.solve_QP(variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate.constraints_jacobian, current_iterate.hessian, d0);
    direction.phase = OPTIMALITY;
    this->number_subproblems_solved++;
    // recompute active set: constraints are active when p-n = 0
    this->recover_l1qp_active_set_(problem, direction, elastic_variables);
    /* remove p and n */
    direction.x.resize(current_iterate.x.size());
    direction.multipliers.lower_bounds.resize(current_iterate.x.size());
    direction.multipliers.upper_bounds.resize(current_iterate.x.size());
    direction.norm = norm_inf(direction.x);

    /* remove contribution of positive part variables */
    for (std::pair<const int, int>& element: elastic_variables.positive) {
        int j = element.first;
        int i = element.second;
        current_iterate.constraints_jacobian[j].erase(i);
        current_iterate.objective_gradient.erase(i);
    }
    /* remove contribution of negative part variables */
    for (std::pair<const int, int>& element: elastic_variables.negative) {
        int j = element.first;
        int i = element.second;
        current_iterate.constraints_jacobian[j].erase(i);
        current_iterate.objective_gradient.erase(i);
    }
    DEBUG << direction;
    return direction;
}

void ActiveSetMethod::recover_l1qp_active_set_(Problem& problem, Direction& direction, const ElasticVariables& elastic_variables) {
    // remove extra variables p and n
    for (unsigned int i = problem.number_variables; i < direction.x.size(); i++) {
        direction.active_set.bounds.at_lower_bound.erase(i);
        direction.active_set.bounds.at_upper_bound.erase(i);
    }
    // constraints: only when p-n = 0
    for (unsigned int j = 0; j < direction.multipliers.constraints.size(); j++) {
        // compute constraint violation
        double constraint_violation = 0.;
        try {
            int i = elastic_variables.positive.at(j);
            constraint_violation += direction.x[i];
        }
        catch (const std::out_of_range& e) {}
        try {
            int i = elastic_variables.negative.at(j);
            constraint_violation += direction.x[i];
        }
        catch (const std::out_of_range& e) {}
        // update active set
        if (0. < constraint_violation) {
            direction.active_set.constraints.at_lower_bound.erase(j);
            direction.active_set.constraints.at_upper_bound.erase(j);
        }
    }
    return;
}

void ActiveSetMethod::generate_elastic_variables_(Problem& problem, ElasticVariables& elastic_variables) {
    // generate elastic variables p and n on the fly to relax the constraints
    int current_index = problem.number_variables;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (-INFINITY < problem.constraint_bounds[j].lb) {
            // nonpositive variable n that captures the negative part of the constraint violation
            elastic_variables.negative[j] = current_index;
            current_index++;
        }
        if (problem.constraint_bounds[j].ub < INFINITY) {
            // nonnegative variable p that captures the positive part of the constraint violation
            elastic_variables.positive[j] = current_index;
            current_index++;
        }
    }
    return;
}

void ActiveSetMethod::compute_l1_linear_objective_(Iterate& current_iterate, ConstraintPartition& constraint_partition) {
    /* objective function: sum of gradients of infeasible constraints */
    SparseGradient objective_gradient;
    for (int j: constraint_partition.infeasible) {
        for (std::pair<int, double> term: current_iterate.constraints_jacobian[j]) {
            int i = term.first;
            double derivative = term.second;

            if (constraint_partition.constraint_feasibility[j] == INFEASIBLE_LOWER) {
                objective_gradient[i] -= derivative;
            }
            else {
                objective_gradient[i] += derivative;
            }
        }
    }
    current_iterate.set_objective_gradient(objective_gradient);
    return;
}

std::vector<double> ActiveSetMethod::generate_l1_multipliers_(Problem& problem, std::vector<double>& current_constraint_multipliers, ConstraintPartition& constraint_partition) {
    std::vector<double> constraint_multipliers(problem.number_constraints);
    for (int j = 0; j < problem.number_constraints; j++) {
        if (constraint_partition.constraint_feasibility[j] == INFEASIBLE_LOWER) {
            constraint_multipliers[j] = 1.;
        }
        else if (constraint_partition.constraint_feasibility[j] == INFEASIBLE_UPPER) {
            constraint_multipliers[j] = -1.;
        }
        else {
            constraint_multipliers[j] = current_constraint_multipliers[j];
        }
    }
    return constraint_multipliers;
}

std::vector<Range> ActiveSetMethod::generate_feasibility_bounds_(Problem& problem, std::vector<double>& current_constraints, ConstraintPartition& constraint_partition) {
    std::vector<Range> constraints_bounds(problem.number_constraints);
    for (int j = 0; j < problem.number_constraints; j++) {
        double lb, ub;
        if (constraint_partition.constraint_feasibility[j] == INFEASIBLE_LOWER) {
            lb = -INFINITY;
            ub = problem.constraint_bounds[j].lb - current_constraints[j];
        }
        else if (constraint_partition.constraint_feasibility[j] == INFEASIBLE_UPPER) {
            lb = problem.constraint_bounds[j].ub - current_constraints[j];
            ub = INFINITY;
        }
        else { // FEASIBLE
            lb = problem.constraint_bounds[j].lb - current_constraints[j];
            ub = problem.constraint_bounds[j].ub - current_constraints[j];
        }
        constraints_bounds[j] = {lb, ub};
    }
    return constraints_bounds;
}

/* LP */

Direction ActiveSetMethod::compute_lp_step_(Problem& problem, QPSolver& solver, Iterate& current_iterate, double trust_region_radius) {
    DEBUG << "Current point: "; print_vector(DEBUG, current_iterate.x);
    DEBUG << "Current constraint multipliers: "; print_vector(DEBUG, current_iterate.multipliers.constraints);
    DEBUG << "Current lb multipliers: "; print_vector(DEBUG, current_iterate.multipliers.lower_bounds);
    DEBUG << "Current ub multipliers: "; print_vector(DEBUG, current_iterate.multipliers.upper_bounds);
    
    /* bounds of the variables */
    std::vector<Range> variables_bounds = this->generate_variables_bounds_(problem, current_iterate, trust_region_radius);

    /* bounds of the linearized constraints */
    std::vector<Range> constraints_bounds = Subproblem::generate_constraints_bounds(problem, current_iterate.constraints);

    /* generate the initial point */
    std::vector<double> d0(current_iterate.x.size()); // = {0.}

    /* solve the QP */
    Direction direction = solver.solve_LP(variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate.constraints_jacobian, d0);
    direction.objective_multiplier = problem.objective_sign;
    direction.phase = OPTIMALITY;
//    direction.predicted_reduction = [&](double step_length) {
//        return this->compute_lp_predicted_reduction_(direction, step_length);
//    };
    direction.predicted_reduction = this->compute_lp_predicted_reduction_;
    this->number_subproblems_solved++;
    DEBUG << direction;
    return direction;
}

double ActiveSetMethod::compute_lp_predicted_reduction_(Problem& /*problem*/, Iterate& /*current_iterate*/, Direction& direction, double step_length) {
    // the predicted reduction is linear in the step length 
    return -step_length*direction.objective;
}

Direction ActiveSetMethod::compute_l1lp_step_(Problem& problem, QPSolver& solver, Iterate& current_iterate, Direction& phase_2_direction, double trust_region_radius) {
    DEBUG << "\nCreating the restoration problem with " << phase_2_direction.constraint_partition.infeasible.size() << " infeasible constraints\n";
    
    /* compute the objective */
    this->compute_l1_linear_objective_(current_iterate, phase_2_direction.constraint_partition);

    /* bounds of the variables */
    std::vector<Range> variables_bounds = this->generate_variables_bounds_(problem, current_iterate, trust_region_radius);

    /* bounds of the linearized constraints */
    std::vector<Range> constraints_bounds = this->generate_feasibility_bounds_(problem, current_iterate.constraints, phase_2_direction.constraint_partition);

    /* generate the initial point */
    std::vector<double> d0 = phase_2_direction.x;

    /* solve the QP */
    Direction direction = solver.solve_LP(variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate.constraints_jacobian, d0);
    direction.objective_multiplier = 0.;
    direction.phase = RESTORATION;
    direction.constraint_partition = phase_2_direction.constraint_partition;
//    direction.predicted_reduction = [&](double step_length) {
//        return this->compute_lp_predicted_reduction_(direction, step_length);
//    };
    direction.predicted_reduction = this->compute_lp_predicted_reduction_;
    this->number_subproblems_solved++;
    DEBUG << direction;
    return direction;
}
