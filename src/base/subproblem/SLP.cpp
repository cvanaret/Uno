#include <cmath>
#include <map>
#include "SLP.hpp"
#include "Constraint.hpp"
#include "Utils.hpp"
#include "Logger.hpp"

SLP::SLP(LPSolver& solver) : Subproblem(), solver(solver) {
}

Iterate SLP::initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers, int number_variables, int number_constraints, bool use_trust_region) {
    Iterate first_iterate(problem, x, multipliers);
    /* compute the optimality and feasibility measures of the initial point */
    this->compute_measures(problem, first_iterate);

    /* allocate the QP solver */
    this->solver.allocate(number_variables, number_constraints);
    return first_iterate;
}

SubproblemSolution SLP::compute_optimality_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds) {
    /* compute first-order information */
    current_iterate.compute_objective_gradient(problem);
    current_iterate.compute_constraints_jacobian(problem);

    /* bounds of the linearized constraints */
    std::vector<Range> constraints_bounds(problem.number_constraints);
    for (int j = 0; j < problem.number_constraints; j++) {
        double lb = problem.constraints_bounds[j].lb - current_iterate.constraints[j];
        double ub = problem.constraints_bounds[j].ub - current_iterate.constraints[j];
        constraints_bounds[j] = {lb, ub};
    }

    /* generate the initial solution */
    std::vector<double> d0(variables_bounds.size()); // = {0.}

    /* solve the QP */
    SubproblemSolution solution = this->solver.solve_LP(variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate.constraints_jacobian, d0);
    solution.phase_1_required = this->phase_1_required(solution);
    this->number_subproblems_solved++;
    return solution;
}

SubproblemSolution SLP::compute_infeasibility_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds, SubproblemSolution& phase_II_solution) {
    DEBUG << "Creating the restoration problem with " << phase_II_solution.constraint_partition.infeasible_set.size() << " infeasible constraints\n";

    /* update the multipliers of the general constraints */
    std::vector<double> constraint_multipliers(problem.number_constraints);
    for (int j = 0; j < problem.number_constraints; j++) {
        if (phase_II_solution.constraint_partition.constraint_status[j] == INFEASIBLE_LOWER) {
            constraint_multipliers[j] = 1.;
        } else if (phase_II_solution.constraint_partition.constraint_status[j] == INFEASIBLE_UPPER) {
            constraint_multipliers[j] = -1.;
        } else {
            constraint_multipliers[j] = current_iterate.multipliers.constraints[j];
        }
    }
    /* compute first-order information */
    current_iterate.compute_objective_gradient(problem);
    current_iterate.compute_constraints_jacobian(problem);

    /* bounds of the linearized constraints */
    std::vector<Range> constraints_bounds(problem.number_constraints);
    for (int j = 0; j < problem.number_constraints; j++) {
        double lb, ub;
        if (phase_II_solution.constraint_partition.constraint_status[j] == INFEASIBLE_LOWER) {
            lb = -INFINITY;
            ub = problem.constraints_bounds[j].lb - current_iterate.constraints[j];
        } else if (phase_II_solution.constraint_partition.constraint_status[j] == INFEASIBLE_UPPER) {
            lb = problem.constraints_bounds[j].ub - current_iterate.constraints[j];
            ub = INFINITY;
        } else { // FEASIBLE
            lb = problem.constraints_bounds[j].lb - current_iterate.constraints[j];
            ub = problem.constraints_bounds[j].ub - current_iterate.constraints[j];
        }
        constraints_bounds[j] = {lb, ub};
    }

    /* compute the objective */
    this->set_infeasibility_objective_(problem, current_iterate, phase_II_solution.constraint_partition);

    /* generate the initial solution */
    std::vector<double> d0 = phase_II_solution.x;

    /* solve the QP */
    SubproblemSolution solution = this->solver.solve_LP(variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate.constraints_jacobian, d0);
    this->number_subproblems_solved++;
    return solution;
}

double SLP::compute_predicted_reduction(Iterate& current_iterate, SubproblemSolution& solution, double step_length) {
    if (step_length == 1.) {
        return -solution.objective;
    } else {
        double linear_term = dot(solution.x, current_iterate.objective_gradient);
        return -step_length*linear_term;
    }
}

/* additional variables */
SubproblemSolution SLP::compute_l1_penalty_step(Problem& problem, Iterate& current_iterate, std::vector<Range>& variables_bounds, double penalty_parameter, PenaltyDimensions penalty_dimensions) {
    /* generate the QP */
    //QP qp = this->generate_l1_penalty_qp_(problem, current_iterate, variables_bounds, penalty_parameter, penalty_dimensions);

    /* generate the initial solution */
    //std::vector<double> d0(problem.number_variables); // = {0.}

    /* solve the QP */
    //SubproblemSolution solution = this->solver.solve(variables_bounds, qp.constraints_bounds, qp.linear_objective, qp.constraints, current_iterate.hessian, d0);


    ActiveSet active_set;
    ConstraintPartition constraint_partition;

    SubproblemSolution solution(current_iterate.x, current_iterate.multipliers, active_set, constraint_partition);
    this->number_subproblems_solved++;
    return solution;
}

/* private methods */

void SLP::set_infeasibility_objective_(Problem& problem, Iterate& current_iterate, ConstraintPartition& constraint_partition) {
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
            } else {
                objective_gradient[variable_index] += derivative;
            }
        }
    }
    current_iterate.set_objective_gradient(objective_gradient);
    return;
}

void SLP::compute_measures(Problem& problem, Iterate& iterate) {
    iterate.feasibility_measure = iterate.residual;
    iterate.optimality_measure = iterate.objective;
    return;
}

bool SLP::phase_1_required(SubproblemSolution& solution) {
    return solution.status == INFEASIBLE;
}