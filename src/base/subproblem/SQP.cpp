#include <cmath>
#include <map>
#include "SQP.hpp"
#include "Constraint.hpp"
#include "Utils.hpp"
#include "Logger.hpp"

SQP::SQP(QPSolver& solver) : Subproblem(), solver(solver) {
}

Iterate SQP::initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers, int number_variables, bool /*use_trust_region*/) {
    // register the original bounds
    this->subproblem_variables_bounds = problem.variables_bounds;

    Iterate first_iterate(problem, x, multipliers);
    /* compute the optimality and feasibility measures of the initial point */
    this->compute_measures(problem, first_iterate);

    /* allocate the QP solver */
    this->solver.allocate(number_variables, problem.number_constraints);

    /* compute least-square multipliers */
    if (0 < problem.number_constraints) {
        first_iterate.compute_constraints_jacobian(problem);
        first_iterate.multipliers.constraints = Subproblem::compute_least_square_multipliers(problem, first_iterate, multipliers.constraints, 1e4);
    }
    return first_iterate;
}

SubproblemSolution SQP::compute_optimality_step(Problem& problem, Iterate& current_iterate, double trust_region_radius) {
    /* evaluate the functions at the current iterate */
    this->evaluate_optimality_iterate(problem, current_iterate);

    /* bounds of the variables */
    std::vector<Range> variables_bounds = this->generate_variables_bounds(current_iterate, trust_region_radius);

    /* bounds of the linearized constraints */
    std::vector<Range> constraints_bounds = Subproblem::generate_constraints_bounds(problem, current_iterate.constraints);

    /* generate the initial point */
    std::vector<double> d0(current_iterate.x.size()); // = {0.}

    /* solve the QP */
    SubproblemSolution solution = this->solve_subproblem(variables_bounds, constraints_bounds, current_iterate, d0);
    solution.phase_1_required = this->phase_1_required(solution);
    this->number_subproblems_solved++;
    return solution;
}

SubproblemSolution SQP::compute_infeasibility_step(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution, double trust_region_radius) {
    DEBUG << "Creating the restoration problem with " << phase_II_solution.constraint_partition.infeasible_set.size() << " infeasible constraints\n";
    
    /* evaluate the functions at the current iterate */
    this->evaluate_feasibility_iterate(problem, current_iterate, phase_II_solution);

    /* bounds of the variables */
    std::vector<Range> variables_bounds = this->generate_variables_bounds(current_iterate, trust_region_radius);

    /* bounds of the linearized constraints */
    std::vector<Range> constraints_bounds = this->generate_feasibility_bounds(problem, current_iterate.constraints, phase_II_solution.constraint_partition);

    /* compute the objective */
    this->set_feasibility_objective_(current_iterate, phase_II_solution.constraint_partition);

    /* generate the initial point */
    std::vector<double> d0 = phase_II_solution.x;

    /* solve the QP */
    SubproblemSolution solution = this->solve_subproblem(variables_bounds, constraints_bounds, current_iterate, d0);
    this->number_subproblems_solved++;
    return solution;
}

double SQP::compute_predicted_reduction(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution, double step_length) {
    if (step_length == 1.) {
        /* full step */
        return -solution.objective;
    }
    else {
        /* the predicted reduction is a quadratic in the step length */
        double linear_term = dot(solution.x, current_iterate.objective_gradient);
        double quadratic_term = current_iterate.hessian.quadratic_product(solution.x, solution.x) / 2.;
        return -step_length * (linear_term + step_length * quadratic_term);
    }
}

void SQP::compute_measures(Problem& problem, Iterate& iterate) {
    iterate.feasibility_measure = iterate.residual;
    iterate.optimality_measure = iterate.objective;
    return;
}

bool SQP::phase_1_required(SubproblemSolution& solution) {
    return (solution.status == INFEASIBLE);
}

/* private methods */

void SQP::evaluate_optimality_iterate(Problem& problem, Iterate& current_iterate) {
    /* compute first- and second-order information */
    current_iterate.compute_hessian(problem, problem.objective_sign, current_iterate.multipliers.constraints);
    current_iterate.compute_objective_gradient(problem);
    current_iterate.compute_constraints_jacobian(problem);
    return;
}

void SQP::evaluate_feasibility_iterate(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution) {
    /* update the multipliers of the general constraints */
    std::vector<double> constraint_multipliers = this->generate_feasibility_multipliers(problem, current_iterate.multipliers.constraints, phase_II_solution.constraint_partition);
    /* compute first- and second-order information */
    double objective_multiplier = 0.;
    current_iterate.compute_hessian(problem, objective_multiplier, constraint_multipliers);
    current_iterate.compute_constraints_jacobian(problem);
    return;
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

void SQP::set_feasibility_objective_(Iterate& current_iterate, ConstraintPartition& constraint_partition) {
    /* objective function: add the gradients of infeasible constraints */
    std::map<int, double> objective_gradient;
    for (int j : constraint_partition.infeasible_set) {
        for (std::pair<int, double> term : current_iterate.constraints_jacobian[j]) {
            int i = term.first;
            double derivative = term.second;

            if (constraint_partition.constraint_status[j] == INFEASIBLE_LOWER) {
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

SubproblemSolution SQP::solve_subproblem(std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds, Iterate& current_iterate, std::vector<double>& d0) {
    return this->solver.solve_QP(variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate.constraints_jacobian, current_iterate.hessian, d0);
}