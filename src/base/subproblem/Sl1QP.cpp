#include <cmath>
#include <map>
#include "Sl1QP.hpp"
#include "Constraint.hpp"
#include "Utils.hpp"
#include "Logger.hpp"

Sl1QP::Sl1QP(QPSolver& solver) : Subproblem(), solver(solver) {
}

Iterate Sl1QP::initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers, int /*number_variables*/, int number_constraints, std::vector<Range>& variables_bounds, bool use_trust_region) {
    // register the original bounds
    this->reformulated_variables_bounds = variables_bounds;
    
    // number of variables in the reformulated problem: |x| + |s| + |p| + |n|
    int number_variables = problem.number_variables + problem.inequality_constraints.size() + 2*problem.number_constraints;
    this->reformulated_variables_bounds.resize(number_variables);
    std::vector<double> reformulated_x(number_variables);
    for (int i = 0; i < problem.number_variables; i++) {
        reformulated_x[i] = x[i];
    }
    Iterate first_iterate(problem, reformulated_x, multipliers);

    // reformulate the problem by introducing slack and relaxation variables
    int current_index = problem.number_variables;
    for (int j = 0; j < problem.number_constraints; j++) {
        if (problem.constraint_status[j] != EQUAL_BOUNDS) {
            // inequality constraints: introduce a slack s
            this->slack_variables[j] = current_index;
            first_iterate.x[current_index] = first_iterate.constraints[j];
            this->reformulated_variables_bounds[current_index] = problem.constraints_bounds[j];
            current_index++;
        }
        // from now on, all constraints are handled as equalities
        // add a nonnegative variable p that captures the positive part of an equality
        this->positive_part_variables[j] = current_index;
        first_iterate.x[current_index] = 0.;
        this->reformulated_variables_bounds[current_index] = {0., INFINITY};
        current_index++;
        // add a nonnegative variable p that captures the positive part of an equality
        this->negative_part_variables[j] = current_index;
        first_iterate.x[current_index] = 0.;
        this->reformulated_variables_bounds[current_index] = {0., INFINITY};
        current_index++;
    }

    /* compute the optimality and feasibility measures of the initial point */
    this->compute_measures(problem, first_iterate);
    
    /* allocate the QP solver */
    this->solver.allocate(number_variables, number_constraints);
    return first_iterate;
}

SubproblemSolution Sl1QP::compute_optimality_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds) {
    /* compute first- and second-order information */
    current_iterate.compute_hessian(problem, problem.objective_sign, current_iterate.multipliers.constraints);
    current_iterate.compute_objective_gradient(problem);
    current_iterate.compute_constraints_jacobian(problem);
    
    /* add contribution of slack variables */
    for (std::pair<const int, int>& element : problem.inequality_constraints) {
        int j = element.first;
        int slack_index = problem.number_variables + element.second;
        first_iterate.constraints_jacobian[j][slack_index] = -1.;
    }
    
    /* bounds of the linearized constraints */
    std::vector<Range> constraints_bounds = Subproblem::generate_constraints_bounds(problem, current_iterate.constraints);

    /* generate the initial point */
    std::vector<double> d0(variables_bounds.size()); // = {0.}

    /* solve the QP */
    SubproblemSolution solution = this->solver.solve_QP(variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate.constraints_jacobian, current_iterate.hessian, d0);
    solution.phase_1_required = this->phase_1_required(solution);
    this->number_subproblems_solved++;
    return solution;
}

SubproblemSolution Sl1QP::compute_infeasibility_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds, SubproblemSolution& phase_II_solution) {
    throw std::out_of_range("Sl1QP.compute_infeasibility_step is not implemented, since l1QP are always feasible");
}

double Sl1QP::compute_predicted_reduction(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution, double step_length) {
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

void Sl1QP::compute_measures(Problem& problem, Iterate& iterate) {
    iterate.feasibility_measure = iterate.residual;
    iterate.optimality_measure = iterate.objective;
    return;
}

bool Sl1QP::phase_1_required(SubproblemSolution& solution) {
    return false;
}

/* private methods */

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
//    int current_index = problem.number_variables;
//    int current_constraint = 0;
//    for (int j = 0; j < problem.number_constraints; j++) {
//        if (problem.constraint_status[j] == EQUAL_BOUNDS) {
//            /* a single constraint with both additional variables */
//            std::map<int, double> gradient(current_iterate.constraints_jacobian[j]);
//            gradient[current_index] = -1.;
//            gradient[current_index + 1] = 1.;
//            qp.constraints[current_constraint] = gradient;
//            /* identical bounds */
//            double lb = problem.constraints_bounds[j].lb - current_iterate.constraints[j];
//            double ub = problem.constraints_bounds[j].ub - current_iterate.constraints[j];
//            qp.constraints_bounds[current_constraint] = {lb, ub};
//            current_index += 2;
//            current_constraint++;
//        }
//        if (problem.constraint_status[j] == BOUNDED_BOTH_SIDES || problem.constraint_status[j] == BOUNDED_LOWER) {
//            /* a single constraint with one additional variable */
//            std::map<int, double> gradient(current_iterate.constraints_jacobian[j]);
//            gradient[current_index] = 1.;
//            qp.constraints[current_constraint] = gradient;
//            /* bounds */
//            double lb = problem.constraints_bounds[j].lb - current_iterate.constraints[j];
//            double ub = INFINITY;
//            qp.constraints_bounds[current_constraint] = {lb, ub};
//            current_index++;
//            current_constraint++;
//        }
//        if (problem.constraint_status[j] == BOUNDED_BOTH_SIDES || problem.constraint_status[j] == BOUNDED_UPPER) {
//            /* a single constraint with one additional variable */
//            std::map<int, double> gradient(current_iterate.constraints_jacobian[j]);
//            gradient[current_index] = -1.;
//            qp.constraints[current_constraint] = gradient;
//            /* bounds */
//            double lb = -INFINITY;
//            double ub = problem.constraints_bounds[j].ub - current_iterate.constraints[j];
//            qp.constraints_bounds[current_constraint] = {lb, ub};
//            current_index++;
//            current_constraint++;
//        }
//    }
//    return qp;
//}
