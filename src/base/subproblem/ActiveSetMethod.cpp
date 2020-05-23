#include <cmath>
#include <map>
#include "ActiveSetMethod.hpp"
#include "Constraint.hpp"
#include "Utils.hpp"
#include "Logger.hpp"

ActiveSetMethod::ActiveSetMethod(Problem& problem, std::shared_ptr<QPSolver> solver, bool scale_residuals):
Subproblem("l1", problem.variables_bounds, scale_residuals), solver(solver) {
}

Iterate ActiveSetMethod::evaluate_initial_point(Problem& problem, std::vector<double>& x, Multipliers& multipliers) {
    Iterate first_iterate(x, multipliers);
    /* compute the optimality and feasibility measures of the initial point */
    this->compute_optimality_measures(problem, first_iterate);
    return first_iterate;
}

std::vector<Range> ActiveSetMethod::generate_variables_bounds(Problem& /*problem*/, Iterate& current_iterate, double trust_region_radius) {
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

void ActiveSetMethod::compute_infeasibility_measures(Problem& problem, Iterate& iterate, SubproblemSolution& solution) {
    iterate.compute_constraints(problem);
    // feasibility measure: residual of all constraints
    iterate.feasibility_measure = problem.compute_constraint_residual(iterate.constraints, this->residual_norm);
    // optimality measure: residual of linearly infeasible constraints
    iterate.optimality_measure = problem.compute_constraint_residual(iterate.constraints, solution.constraint_partition.infeasible, this->residual_norm);
    return;
}

/* QP */

SubproblemSolution ActiveSetMethod::compute_qp_step(Problem& problem, Iterate& current_iterate, double trust_region_radius) {
    DEBUG << "Current point: "; print_vector(DEBUG, current_iterate.x);
    DEBUG << "Current constraint multipliers: "; print_vector(DEBUG, current_iterate.multipliers.constraints);
    DEBUG << "Current lb multipliers: "; print_vector(DEBUG, current_iterate.multipliers.lower_bounds);
    DEBUG << "Current ub multipliers: "; print_vector(DEBUG, current_iterate.multipliers.upper_bounds);
    
    /* bounds of the variables */
    std::vector<Range> variables_bounds = this->generate_variables_bounds(problem, current_iterate, trust_region_radius);

    /* bounds of the linearized constraints */
    std::vector<Range> constraints_bounds = Subproblem::generate_constraints_bounds(problem, current_iterate.constraints);

    /* generate the initial point */
    std::vector<double> d0(variables_bounds.size()); // = {0.}

    /* solve the QP */
    SubproblemSolution solution = this->solver->solve_QP(variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate.constraints_jacobian, current_iterate.hessian, d0);
    solution.phase = OPTIMALITY;
    this->number_subproblems_solved++;
    DEBUG << solution;
    return solution;
}

double ActiveSetMethod::compute_qp_predicted_reduction(Iterate& current_iterate, SubproblemSolution& solution, double step_length) {
    // the predicted reduction is quadratic in the step length
    if (step_length == 1.) {
        return -solution.objective;
    }
    else {
        double linear_term = dot(solution.x, current_iterate.objective_gradient);
        double quadratic_term = current_iterate.hessian.quadratic_product(solution.x, solution.x) / 2.;
        return -step_length * (linear_term + step_length * quadratic_term);
    }
}

SubproblemSolution ActiveSetMethod::compute_feasibility_qp_step(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution, double trust_region_radius) {
    DEBUG << "\nCreating the restoration problem with " << phase_II_solution.constraint_partition.infeasible.size() << " infeasible constraints\n";

    /* compute the objective */
    this->compute_linear_feasibility_objective(current_iterate, phase_II_solution.constraint_partition);

    /* bounds of the variables */
    std::vector<Range> variables_bounds = this->generate_variables_bounds(problem, current_iterate, trust_region_radius);

    /* bounds of the linearized constraints */
    std::vector<Range> constraints_bounds = this->generate_feasibility_bounds(problem, current_iterate.constraints, phase_II_solution.constraint_partition);

    /* generate the initial point */
    std::vector<double> d0 = phase_II_solution.x;

    /* solve the QP */
    SubproblemSolution solution = this->solver->solve_QP(variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate.constraints_jacobian, current_iterate.hessian, d0);
    solution.objective_multiplier = 0.;
    solution.phase = RESTORATION;
    solution.constraint_partition = phase_II_solution.constraint_partition;
    this->number_subproblems_solved++;
    DEBUG << solution;
    return solution;
}

void ActiveSetMethod::compute_linear_feasibility_objective(Iterate& current_iterate, ConstraintPartition& constraint_partition) {
    /* objective function: sum of gradients of infeasible constraints */
    std::map<int, double> objective_gradient;
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

std::vector<double> ActiveSetMethod::generate_feasibility_multipliers(Problem& problem, std::vector<double>& current_constraint_multipliers, ConstraintPartition& constraint_partition) {
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

std::vector<Range> ActiveSetMethod::generate_feasibility_bounds(Problem& problem, std::vector<double>& current_constraints, ConstraintPartition& constraint_partition) {
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

SubproblemSolution ActiveSetMethod::compute_lp_step(Problem& problem, Iterate& current_iterate, double trust_region_radius) {
    DEBUG << "Current point: "; print_vector(DEBUG, current_iterate.x);
    DEBUG << "Current constraint multipliers: "; print_vector(DEBUG, current_iterate.multipliers.constraints);
    DEBUG << "Current lb multipliers: "; print_vector(DEBUG, current_iterate.multipliers.lower_bounds);
    DEBUG << "Current ub multipliers: "; print_vector(DEBUG, current_iterate.multipliers.upper_bounds);
    
    /* bounds of the variables */
    std::vector<Range> variables_bounds = this->generate_variables_bounds(problem, current_iterate, trust_region_radius);

    /* bounds of the linearized constraints */
    std::vector<Range> constraints_bounds = Subproblem::generate_constraints_bounds(problem, current_iterate.constraints);

    /* generate the initial point */
    std::vector<double> d0(current_iterate.x.size()); // = {0.}

    /* solve the QP */
    SubproblemSolution solution = this->solver->solve_LP(variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate.constraints_jacobian, d0);
    solution.objective_multiplier = problem.objective_sign;
    solution.phase = OPTIMALITY;
    solution.predicted_reduction = [&](double step_length) {
        return this->compute_lp_predicted_reduction(solution, step_length);
    };
    this->number_subproblems_solved++;
    DEBUG << solution;
    return solution;
}

double ActiveSetMethod::compute_lp_predicted_reduction(SubproblemSolution& solution, double step_length) {
    // the predicted reduction is linear in the step length 
    return -step_length*solution.objective;
}

SubproblemSolution ActiveSetMethod::compute_feasibility_lp_step(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution, double trust_region_radius) {
    DEBUG << "\nCreating the restoration problem with " << phase_II_solution.constraint_partition.infeasible.size() << " infeasible constraints\n";
    
    /* compute the objective */
    this->compute_linear_feasibility_objective(current_iterate, phase_II_solution.constraint_partition);

    /* bounds of the variables */
    std::vector<Range> variables_bounds = this->generate_variables_bounds(problem, current_iterate, trust_region_radius);

    /* bounds of the linearized constraints */
    std::vector<Range> constraints_bounds = this->generate_feasibility_bounds(problem, current_iterate.constraints, phase_II_solution.constraint_partition);

    /* generate the initial point */
    std::vector<double> d0 = phase_II_solution.x;

    /* solve the QP */
    SubproblemSolution solution = this->solver->solve_LP(variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate.constraints_jacobian, d0);
    solution.objective_multiplier = 0.;
    solution.phase = RESTORATION;
    solution.constraint_partition = phase_II_solution.constraint_partition;
    solution.predicted_reduction = [&](double step_length) {
        return this->compute_lp_predicted_reduction(solution, step_length);
    };
    this->number_subproblems_solved++;
    DEBUG << solution;
    return solution;
}