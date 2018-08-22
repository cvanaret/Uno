#include <cmath>
#include <map>
#include "QPApproximation.hpp"
#include "Constraint.hpp"
#include "Utils.hpp"
#include "Logger.hpp"

QPApproximation::QPApproximation(QPSolver& solver) : Subproblem(), solver(solver) {
}

Iterate QPApproximation::initialize(Problem& problem, std::vector<double>& x, std::vector<double>& bound_multipliers, std::vector<double>& constraint_multipliers, int number_variables, int number_constraints, bool use_trust_region) {
    Iterate first_iterate(problem, x, bound_multipliers, constraint_multipliers);
    /* compute the optimality and feasibility measures of the initial point */
    this->compute_measures(problem, first_iterate);

    /* allocate the QP solver */
    this->solver.allocate(number_variables, number_constraints);
    return first_iterate;
}

LocalSolution QPApproximation::compute_optimality_step(Problem& problem, Iterate& current_iterate, double radius) {
    /* generate the QP */
    QP qp = this->generate_optimality_qp_(problem, current_iterate, radius);
    DEBUG << qp;

    /* generate the initial solution */
    std::vector<double> d0(qp.number_variables); // = {0.}

    /* solve the QP */
    LocalSolution solution = this->solver.solve(qp, d0);
    this->number_subproblems_solved++;

    /* keep multipliers delta */
    for (int j = 0; j < problem.number_constraints; j++) {
        solution.constraint_multipliers[j] -= current_iterate.constraint_multipliers[j];
    }

    double linear_term = dot(solution.x, qp.objective);
    double quadratic_term = current_iterate.hessian.quadratic_product(solution.x, solution.x) / 2.;
    solution.objective_terms = {linear_term, quadratic_term};

    return solution;
}

LocalSolution QPApproximation::compute_infeasibility_step(Problem& problem, Iterate& current_iterate, double radius, LocalSolution& phase_II_solution) {
    /* generate the QP */
    QP qp = this->generate_infeasibility_qp_(problem, current_iterate, radius, phase_II_solution.constraint_partition);
    DEBUG << qp;

    /* generate the initial solution */
    std::vector<double> d0 = phase_II_solution.x;

    /* solve the QP */
    LocalSolution solution = this->solver.solve(qp, d0);
    this->number_subproblems_solved++;

    /* keep multipliers delta */
    for (int j = 0; j < problem.number_constraints; j++) {
        solution.constraint_multipliers[j] -= current_iterate.constraint_multipliers[j];
    }

    double linear_term = dot(solution.x, qp.objective);
    double quadratic_term = current_iterate.hessian.quadratic_product(solution.x, solution.x) / 2.;
    solution.objective_terms = {linear_term, quadratic_term};

    return solution;
}

/* additional variables */
LocalSolution QPApproximation::compute_l1_penalty_step(Problem& problem, Iterate& current_iterate, double radius, double penalty_parameter, PenaltyDimensions penalty_dimensions) {
    /* generate the QP */
    QP qp = this->generate_l1_penalty_qp_(problem, current_iterate, radius, penalty_parameter, penalty_dimensions);
    DEBUG << qp;

    /* generate the initial solution */
    std::vector<double> d0(qp.number_variables); // = {0.}

    /* solve the QP */
    LocalSolution solution = this->solver.solve(qp, d0);
    this->number_subproblems_solved++;

    /* keep multipliers delta */
    for (int j = 0; j < problem.number_constraints; j++) {
        solution.constraint_multipliers[j] -= current_iterate.constraint_multipliers[j];
    }

    double linear_term = dot(solution.x, qp.objective);
    double quadratic_term = current_iterate.hessian.quadratic_product(solution.x, solution.x) / 2.;
    solution.objective_terms = {linear_term, quadratic_term};

    return solution;
}

QP QPApproximation::generate_qp_(Problem& problem, Iterate& current_iterate, double radius) {
    /* initialize the QP */
    QP qp(problem.number_variables, problem.number_constraints, current_iterate.hessian);

    /* bound constraints intersected with trust region  */
    for (int i = 0; i < problem.number_variables; i++) {
        qp.variable_lb[i] = std::max(-radius, problem.variable_lb[i] - current_iterate.x[i]);
        qp.variable_ub[i] = std::min(radius, problem.variable_ub[i] - current_iterate.x[i]);
    }

    /* compute the constraints */
    this->set_constraints_(problem, qp, current_iterate);

    return qp;
}

QP QPApproximation::generate_optimality_qp_(Problem& problem, Iterate& current_iterate, double radius) {
    DEBUG << "Creating the optimality problem\n";

    /* compute the Lagrangian Hessian */
    current_iterate.compute_hessian(problem, problem.objective_sign, current_iterate.constraint_multipliers);

    /* initialize the QP */
    QP qp = this->generate_qp_(problem, current_iterate, radius);

    /* bounds of the linearized constraints */
    for (int j = 0; j < qp.number_constraints; j++) {
        qp.constraint_lb[j] = problem.constraint_lb[j] - current_iterate.constraints[j];
        qp.constraint_ub[j] = problem.constraint_ub[j] - current_iterate.constraints[j];
    }

    /* compute the objective */
    this->set_optimality_objective_(problem, qp, current_iterate);

    return qp;
}

QP QPApproximation::generate_infeasibility_qp_(Problem& problem, Iterate& current_iterate, double radius, ConstraintPartition& constraint_partition) {
    int number_infeasible = constraint_partition.infeasible_set.size();
    DEBUG << "Creating the restoration problem with " << number_infeasible << " infeasible constraints\n";

    /* update the multipliers of the general constraints */
    std::vector<double> constraint_multipliers(problem.number_constraints);
    for (int j = 0; j < problem.number_constraints; j++) {
        if (constraint_partition.constraint_status[j] == INFEASIBLE_LOWER) {
            constraint_multipliers[j] = 1.;
        }
        else if (constraint_partition.constraint_status[j] == INFEASIBLE_UPPER) {
            constraint_multipliers[j] = -1.;
        }
        else {
            constraint_multipliers[j] = current_iterate.constraint_multipliers[j];
        }
    }
    /* compute the Lagrangian Hessian */
    double objective_multiplier = 0.;
    current_iterate.compute_hessian(problem, objective_multiplier, constraint_multipliers);

    /* initialize the QP */
    QP qp = this->generate_qp_(problem, current_iterate, radius);

    /* bounds of the linearized constraints */
    for (int j = 0; j < qp.number_constraints; j++) {
        if (constraint_partition.constraint_status[j] == INFEASIBLE_LOWER) {
            qp.constraint_lb[j] = -INFINITY;
            qp.constraint_ub[j] = problem.constraint_lb[j] - current_iterate.constraints[j];
        }
        else if (constraint_partition.constraint_status[j] == INFEASIBLE_UPPER) {
            qp.constraint_lb[j] = problem.constraint_ub[j] - current_iterate.constraints[j];
            qp.constraint_ub[j] = INFINITY;
        }
        else { // FEASIBLE
            qp.constraint_lb[j] = problem.constraint_lb[j] - current_iterate.constraints[j];
            qp.constraint_ub[j] = problem.constraint_ub[j] - current_iterate.constraints[j];
        }
    }

    /* compute the objective */
    this->set_infeasibility_objective_(problem, qp, current_iterate, constraint_partition);
    return qp;
}

QP QPApproximation::generate_l1_penalty_qp_(Problem& problem, Iterate& current_iterate, double radius, double penalty_parameter, PenaltyDimensions penalty_dimensions) {
    int number_variables = problem.number_variables + penalty_dimensions.number_additional_variables;
    int number_constraints = penalty_dimensions.number_constraints;

    /* compute the Lagrangian Hessian from scratch */
    current_iterate.is_hessian_computed = false;
    double objective_multiplier = penalty_parameter;
    current_iterate.compute_hessian(problem, objective_multiplier, current_iterate.constraint_multipliers);

    /* initialize the QP */
    QP qp(number_variables, number_constraints, current_iterate.hessian);

    /* bounds of original variables intersected with trust region  */
    for (int i = 0; i < problem.number_variables; i++) {
        qp.variable_lb[i] = std::max(-radius, problem.variable_lb[i] - current_iterate.x[i]);
        qp.variable_ub[i] = std::min(radius, problem.variable_ub[i] - current_iterate.x[i]);
    }
    /* bounds of additional variables */
    for (int k = 0; k < penalty_dimensions.number_additional_variables; k++) {
        qp.variable_lb[problem.number_variables + k] = 0.;
        qp.variable_ub[problem.number_variables + k] = INFINITY;
    }

    /* apply the nonzero penalty parameter on the initial objective */
    if (penalty_parameter != 0.) {
        if (!current_iterate.is_objective_gradient_computed) {
            std::map<int, double> objective_gradient = problem.objective_sparse_gradient(current_iterate.x);
            current_iterate.set_objective_gradient(objective_gradient);
        }
        qp.objective = current_iterate.objective_gradient;
        for (std::pair<int, double> term : qp.objective) {
            int index = term.first;
            qp.objective[index] *= penalty_parameter;
        }
    }
    /* add additional variables to the objective */
    for (int k = 0; k < penalty_dimensions.number_additional_variables; k++) {
        qp.objective[problem.number_variables + k] = 1.;
    }

    /* compute the original constraint gradients */
    current_iterate.compute_constraints_jacobian(problem);

    /* add the constraints */
    int current_additional_variable = problem.number_variables;
    int current_constraint = 0;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] == EQUAL_BOUNDS) {
            /* a single constraint with both additional variables */
            std::map<int, double> gradient(current_iterate.constraints_jacobian[j]);
            gradient[current_additional_variable] = -1.;
            gradient[current_additional_variable + 1] = 1.;
            qp.constraints[current_constraint] = gradient;
            /* identical bounds */
            qp.constraint_lb[current_constraint] = problem.constraint_lb[j] - current_iterate.constraints[j];
            qp.constraint_ub[current_constraint] = problem.constraint_ub[j] - current_iterate.constraints[j];
            current_additional_variable += 2;
            current_constraint++;
        }
        if (problem.constraint_status[j] == BOUNDED_BOTH_SIDES || problem.constraint_status[j] == BOUNDED_LOWER) {
            /* a single constraint with one additional variable */
            std::map<int, double> gradient(current_iterate.constraints_jacobian[j]);
            gradient[current_additional_variable] = 1.;
            qp.constraints[current_constraint] = gradient;
            /* bounds */
            qp.constraint_lb[current_constraint] = problem.constraint_lb[j] - current_iterate.constraints[j];
            qp.constraint_ub[current_constraint] = INFINITY;
            current_additional_variable++;
            current_constraint++;
        }
        if (problem.constraint_status[j] == BOUNDED_BOTH_SIDES || problem.constraint_status[j] == BOUNDED_UPPER) {
            /* a single constraint with one additional variable */
            std::map<int, double> gradient(current_iterate.constraints_jacobian[j]);
            gradient[current_additional_variable] = -1.;
            qp.constraints[current_constraint] = gradient;
            /* bounds */
            qp.constraint_lb[current_constraint] = -INFINITY;
            qp.constraint_ub[current_constraint] = problem.constraint_ub[j] - current_iterate.constraints[j];
            current_additional_variable++;
            current_constraint++;
        }
    }
    return qp;
}

void QPApproximation::set_constraints_(Problem& problem, QP& qp, Iterate& current_iterate) {
    /* compute the constraint Jacobian */
    current_iterate.compute_constraints_jacobian(problem);
    qp.constraints = current_iterate.constraints_jacobian;
    return;
}

void QPApproximation::set_optimality_objective_(Problem& problem, QP& qp, Iterate& current_iterate) {
    /* compute the objective Jacobian */
    if (!current_iterate.is_objective_gradient_computed) {
        std::map<int, double> objective_gradient = problem.objective_sparse_gradient(current_iterate.x);
        current_iterate.set_objective_gradient(objective_gradient);
    }
    qp.objective = current_iterate.objective_gradient;
    return;
}

void QPApproximation::set_infeasibility_objective_(Problem& problem, QP& qp, Iterate& current_iterate, ConstraintPartition& constraint_partition) {
    /* objective function: add the gradients of infeasible constraints */
    std::map<int, double> objective_gradient;

    for (unsigned int k = 0; k < constraint_partition.infeasible_set.size(); k++) {
        int j = constraint_partition.infeasible_set[k];
        /* combine into objective_gradient */

        for (std::pair<int, double> term : current_iterate.constraints_jacobian[j]) {
            int variable_index = term.first;
            double derivative = term.second;

            if (constraint_partition.constraint_status[j] == INFEASIBLE_LOWER) {
                objective_gradient[variable_index] -= derivative;
            }
            else {
                objective_gradient[variable_index] += derivative;
            }
        }
    }
    current_iterate.set_objective_gradient(objective_gradient);
    qp.objective = current_iterate.objective_gradient;
    return;
}

void QPApproximation::compute_measures(Problem& problem, Iterate& iterate) {
    iterate.feasibility_measure = iterate.residual;
    iterate.optimality_measure = iterate.objective;
    return;
}