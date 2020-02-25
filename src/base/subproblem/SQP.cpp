#include <cmath>
#include <map>
#include "SQP.hpp"
#include "Constraint.hpp"
#include "Utils.hpp"
#include "Logger.hpp"

SQP::SQP(QPSolver& solver) : Subproblem(), solver(solver) {
}

Iterate SQP::initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers, int number_variables, int number_constraints, bool use_trust_region) {
    Iterate first_iterate(problem, x, multipliers);
    /* compute the optimality and feasibility measures of the initial point */
    this->compute_measures(problem, first_iterate);

    /* allocate the QP solver */
    this->solver.allocate(number_variables, number_constraints);
    return first_iterate;
}

SubproblemSolution SQP::compute_optimality_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds) {
    /* compute first- and second-order information */
    current_iterate.compute_hessian(problem, problem.objective_sign, current_iterate.multipliers.constraints);
    current_iterate.compute_objective_gradient(problem);
    current_iterate.compute_constraints_jacobian(problem);
    
    /* bounds of the linearized constraints */
    std::vector<Range> constraints_bounds = this->generate_optimality_bounds(problem, current_iterate.constraints);

    /* generate the initial solution */
    std::vector<double> d0(variables_bounds.size()); // = {0.}

    /* solve the QP */
    SubproblemSolution solution = this->solver.solve_QP(variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate.constraints_jacobian, current_iterate.hessian, d0);
    solution.phase_1_required = this->phase_1_required(solution);
    this->number_subproblems_solved++;
    return solution;
}

SubproblemSolution SQP::compute_infeasibility_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds, SubproblemSolution& phase_II_solution) {
    /* update the multipliers of the general constraints */
    std::vector<double> constraint_multipliers = this->generate_feasibility_multipliers(problem, current_iterate.multipliers.constraints, phase_II_solution.constraint_partition);
    
    /* compute first- and second-order information */
    double objective_multiplier = 0.;
    current_iterate.compute_hessian(problem, objective_multiplier, constraint_multipliers);
    current_iterate.compute_objective_gradient(problem);
    current_iterate.compute_constraints_jacobian(problem);

    /* bounds of the linearized constraints */
    std::vector<Range> constraints_bounds = this->generate_feasibility_bounds(problem, current_iterate.constraints, phase_II_solution.constraint_partition);

    /* compute the objective */
    this->set_feasibility_objective_(problem, current_iterate, phase_II_solution.constraint_partition);

    /* generate the initial solution */
    std::vector<double> d0 = phase_II_solution.x;

    /* solve the QP */
    SubproblemSolution solution = this->solver.solve_QP(variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate.constraints_jacobian, current_iterate.hessian, d0);
    this->number_subproblems_solved++;
    return solution;
}

double SQP::compute_predicted_reduction(Iterate& current_iterate, SubproblemSolution& solution, double step_length) {
    if (step_length == 1.) {
        /* full step */
        return -solution.objective;
    }
    else {
        /* the predicted reduction is a quadratic in the step length */
        double linear_term = dot(solution.x, current_iterate.objective_gradient);
        double quadratic_term = current_iterate.hessian.quadratic_product(solution.x, solution.x) / 2.;
        return -step_length*(linear_term + step_length * quadratic_term);
    } 
}

/* additional variables */
SubproblemSolution SQP::compute_l1_penalty_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds, double penalty_parameter, PenaltyDimensions penalty_dimensions) {
    /* generate the initial solution */
    std::vector<double> d0(problem.number_variables); // = {0.}

    /* solve the QP */
    // TODO constraints_bounds
    SubproblemSolution solution = this->solver.solve_QP(variables_bounds, variables_bounds, current_iterate.objective_gradient, current_iterate.constraints_jacobian, current_iterate.hessian, d0);
    this->number_subproblems_solved++;
    return solution;
}

/* private methods */

std::vector<Range> SQP::generate_optimality_bounds(Problem& problem, std::vector<double>& current_constraints) {
    std::vector<Range> constraints_bounds(problem.number_constraints);
    for (int j = 0; j < problem.number_constraints; j++) {
        double lb = problem.constraints_bounds[j].lb - current_constraints[j];
        double ub = problem.constraints_bounds[j].ub - current_constraints[j];
        constraints_bounds[j] = {lb, ub};
    }
    return constraints_bounds;
}

std::vector<double> SQP::generate_feasibility_multipliers(Problem& problem, std::vector<double>& current_constraint_multipliers, ConstraintPartition& constraint_partition) {
    std::vector<double> constraint_multipliers(problem.number_constraints);
    for (int j = 0; j < problem.number_constraints; j++) {
        if (constraint_partition.constraint_status[j] == INFEASIBLE_LOWER) {
            constraint_multipliers[j] = 1.;
        }
        else if (constraint_partition.constraint_status[j] == INFEASIBLE_UPPER) {
            constraint_multipliers[j] = -1.;
        }
        else {
            constraint_multipliers[j] = current_constraint_multipliers[j];
        }
    }
    return constraint_multipliers;
}

std::vector<Range> SQP::generate_feasibility_bounds(Problem& problem, std::vector<double>& current_constraints, ConstraintPartition& constraint_partition) {
    std::vector<Range> constraints_bounds(problem.number_constraints);
    for (int j = 0; j < problem.number_constraints; j++) {
        double lb, ub;
        if (constraint_partition.constraint_status[j] == INFEASIBLE_LOWER) {
            lb = -INFINITY;
            ub = problem.constraints_bounds[j].lb - current_constraints[j];
        }
        else if (constraint_partition.constraint_status[j] == INFEASIBLE_UPPER) {
            lb = problem.constraints_bounds[j].ub - current_constraints[j];
            ub = INFINITY;
        }
        else { // FEASIBLE
            lb = problem.constraints_bounds[j].lb - current_constraints[j];
            ub = problem.constraints_bounds[j].ub - current_constraints[j];
        }
        constraints_bounds[j] = {lb, ub};
    }
    return constraints_bounds;
}

void SQP::set_feasibility_objective_(Problem& problem, Iterate& current_iterate, ConstraintPartition& constraint_partition) {
    /* objective function: add the gradients of infeasible constraints */
    std::map<int, double> objective_gradient;
    for (int j: constraint_partition.infeasible_set) {
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
    return;
}

void SQP::compute_measures(Problem& problem, Iterate& iterate) {
    iterate.feasibility_measure = iterate.residual;
    iterate.optimality_measure = iterate.objective;
    return;
}

bool SQP::phase_1_required(SubproblemSolution& solution) {
    return solution.status == INFEASIBLE;
}


// TO MOVE

//QP QPApproximation::generate_l1_penalty_qp_(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds, double penalty_parameter, PenaltyDimensions penalty_dimensions) {
//    int number_variables = problem.number_variables + penalty_dimensions.number_additional_variables;
//    int number_constraints = penalty_dimensions.number_constraints;
//
//    /* compute the Lagrangian Hessian from scratch */
//    current_iterate.is_hessian_computed = false;
//    double objective_multiplier = penalty_parameter;
//    current_iterate.compute_hessian(problem, objective_multiplier, current_iterate.multipliers.constraints);
//
//    /* initialize the QP */
//    QP qp(number_variables, number_constraints);
//    qp.variables_bounds = variables_bounds;
//
//    /* bounds of additional variables */
//    for (int k = 0; k < penalty_dimensions.number_additional_variables; k++) {
//        qp.variables_bounds[problem.number_variables + k] = {0., (double) INFINITY};
//    }
//
//    /* apply the nonzero penalty parameter on the initial objective */
//    if (penalty_parameter != 0.) {
//        if (!current_iterate.is_objective_gradient_computed) {
//            std::map<int, double> objective_gradient = problem.objective_sparse_gradient(current_iterate.x);
//            current_iterate.set_objective_gradient(objective_gradient);
//        }
//        qp.linear_objective = current_iterate.objective_gradient;
//        for (std::pair<int, double> term : qp.linear_objective) {
//            int index = term.first;
//            qp.linear_objective[index] *= penalty_parameter;
//        }
//    }
//    /* add additional variables to the objective */
//    for (int k = 0; k < penalty_dimensions.number_additional_variables; k++) {
//        qp.linear_objective[problem.number_variables + k] = 1.;
//    }
//
//    /* compute the original constraint gradients */
//    current_iterate.compute_constraints_jacobian(problem);
//
//    /* add the constraints */
//    int current_additional_variable = problem.number_variables;
//    int current_constraint = 0;
//    for (int j = 0; j < problem.number_constraints; j++) {
//        if (problem.constraint_status[j] == EQUAL_BOUNDS) {
//            /* a single constraint with both additional variables */
//            std::map<int, double> gradient(current_iterate.constraints_jacobian[j]);
//            gradient[current_additional_variable] = -1.;
//            gradient[current_additional_variable + 1] = 1.;
//            qp.constraints[current_constraint] = gradient;
//            /* identical bounds */
//            double lb = problem.constraints_bounds[j].lb - current_iterate.constraints[j];
//            double ub = problem.constraints_bounds[j].ub - current_iterate.constraints[j];
//            qp.constraints_bounds[current_constraint] = {lb, ub};
//            current_additional_variable += 2;
//            current_constraint++;
//        }
//        if (problem.constraint_status[j] == BOUNDED_BOTH_SIDES || problem.constraint_status[j] == BOUNDED_LOWER) {
//            /* a single constraint with one additional variable */
//            std::map<int, double> gradient(current_iterate.constraints_jacobian[j]);
//            gradient[current_additional_variable] = 1.;
//            qp.constraints[current_constraint] = gradient;
//            /* bounds */
//            double lb = problem.constraints_bounds[j].lb - current_iterate.constraints[j];
//            double ub = INFINITY;
//            qp.constraints_bounds[current_constraint] = {lb, ub};
//            current_additional_variable++;
//            current_constraint++;
//        }
//        if (problem.constraint_status[j] == BOUNDED_BOTH_SIDES || problem.constraint_status[j] == BOUNDED_UPPER) {
//            /* a single constraint with one additional variable */
//            std::map<int, double> gradient(current_iterate.constraints_jacobian[j]);
//            gradient[current_additional_variable] = -1.;
//            qp.constraints[current_constraint] = gradient;
//            /* bounds */
//            double lb = -INFINITY;
//            double ub = problem.constraints_bounds[j].ub - current_iterate.constraints[j];
//            qp.constraints_bounds[current_constraint] = {lb, ub};
//            current_additional_variable++;
//            current_constraint++;
//        }
//    }
//    return qp;
//}