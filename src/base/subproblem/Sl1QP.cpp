#include <cmath>
#include <map>
#include "Sl1QP.hpp"
#include "Constraint.hpp"
#include "Utils.hpp"

Sl1QP::Sl1QP(QPSolver& solver) : Subproblem(), solver(solver) {
}

Iterate Sl1QP::initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers, int /*number_variables*/, bool use_trust_region) {
    // register the original bounds
    this->subproblem_variables_bounds = problem.variables_bounds;
    
    // number of variables in the reformulated problem: |x| + |s| + |p| + |n|
    int number_variables = problem.number_variables + problem.inequality_constraints.size() + 2*problem.number_constraints;
    this->subproblem_variables_bounds.resize(number_variables);
    
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
            first_iterate.x[current_index] = Subproblem::project_variable_in_bounds(first_iterate.constraints[j], problem.constraints_bounds[j]);
            this->subproblem_variables_bounds[current_index] = problem.constraints_bounds[j];
            current_index++;
        }
        // from now on, all constraints are handled as equalities
        // add a nonnegative variable p that captures the positive part of an equality
        this->positive_part_variables[j] = current_index;
        first_iterate.x[current_index] = 0.;
        this->subproblem_variables_bounds[current_index] = {0., INFINITY};
        current_index++;
        // add a nonnegative variable p that captures the positive part of an equality
        this->negative_part_variables[j] = current_index;
        first_iterate.x[current_index] = 0.;
        this->subproblem_variables_bounds[current_index] = {0., INFINITY};
        current_index++;
    }

    /* compute the optimality and feasibility measures of the initial point */
    this->compute_measures(problem, first_iterate);
    
    /* allocate the QP solver */
    this->solver.allocate(number_variables, problem.number_constraints);
    return first_iterate;
}

SubproblemSolution Sl1QP::compute_optimality_step(Problem& problem, Iterate& current_iterate, double trust_region_radius) {
    double penalty_parameter = 1;
    
    /* compute first- and second-order information */
    current_iterate.compute_hessian(problem, problem.objective_sign, current_iterate.multipliers.constraints);
    current_iterate.compute_objective_gradient(problem);
    current_iterate.compute_constraints_jacobian(problem);
    
    /* add contribution of slack variables */
    for (std::pair<const int, int>& element : this->slack_variables) {
        int j = element.first;
        int slack_index = element.second;
        current_iterate.constraints_jacobian[j][slack_index] = -1.;
    }
    /* add contribution of positive part variables */
    for (std::pair<const int, int>& element : this->positive_part_variables) {
        int j = element.first;
        int variable_index = element.second;
        current_iterate.constraints_jacobian[j][variable_index] = -1.;
        current_iterate.objective_gradient[variable_index] = penalty_parameter;
    }
    /* add contribution of negative part variables */
    for (std::pair<const int, int>& element : this->negative_part_variables) {
        int j = element.first;
        int variable_index = element.second;
        current_iterate.constraints_jacobian[j][variable_index] = 1.;
        current_iterate.objective_gradient[variable_index] = penalty_parameter;
    }
    
    /* bounds of the variables */
    std::vector<Range> variables_bounds = this->generate_variables_bounds(problem, current_iterate, trust_region_radius);
    
    /* bounds of the linearized equality constraints */
    std::vector<Range> constraints_bounds(problem.number_constraints);
    for (std::pair<const int, int>& element : this->slack_variables) {
        int j = element.first;
        int slack_index = element.second;
        double bound = -(current_iterate.constraints[j] - current_iterate.x[slack_index]);
        constraints_bounds[j] = {bound, bound};
    }
    for (std::pair<const int, int>& element : problem.equality_constraints) {
        int j = element.first;
        double bound = -(current_iterate.constraints[j] - problem.constraints_bounds[j].lb);
        constraints_bounds[j] = {bound, bound};
    }
    
    DEBUG << "gradient obj: "; print_vector(DEBUG, current_iterate.objective_gradient);
    for (int j = 0; j < problem.number_constraints; j++) {
        DEBUG << "gradient c" << j << ": "; print_vector(DEBUG, current_iterate.constraints_jacobian[j]);
    }
    for (unsigned int i = 0; i < current_iterate.x.size(); i++) {
        DEBUG << "x" << i << " in [" << this->subproblem_variables_bounds[i].lb << ", " << this->subproblem_variables_bounds[i].ub << "]\n";
    }
    for (unsigned int i = 0; i < current_iterate.x.size(); i++) {
        DEBUG << "delta x" << i << " in [" << variables_bounds[i].lb << ", " << variables_bounds[i].ub << "]\n";
    }
    for (int j = 0; j < problem.number_constraints; j++) {
        DEBUG << "c" << j << " in [" << constraints_bounds[j].lb << ", " << constraints_bounds[j].ub << "]\n";
    }

    /* generate the initial point */
    std::vector<double> d0(variables_bounds.size()); // = {0.}

    /* solve the QP */
    SubproblemSolution solution = this->solver.solve_QP(variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate.constraints_jacobian, current_iterate.hessian, d0);
    // recompute active set: constraints are active when p-n = 0
    this->recover_active_set(solution, variables_bounds, constraints_bounds);
    
    solution.phase_1_required = this->phase_1_required(solution);
    this->number_subproblems_solved++;
    return solution;
}

SubproblemSolution Sl1QP::compute_infeasibility_step(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution, double trust_region_radius) {
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

std::vector<Range> Sl1QP::generate_variables_bounds(Problem& problem, Iterate& current_iterate, double trust_region_radius) {
    std::vector<Range> variables_bounds(current_iterate.x.size());
    /* original bounds intersected with trust region  */
    for (int i = 0; i < problem.number_variables; i++) {
        double lb = std::max(-trust_region_radius, this->subproblem_variables_bounds[i].lb - current_iterate.x[i]);
        double ub = std::min(trust_region_radius, this->subproblem_variables_bounds[i].ub - current_iterate.x[i]);
        variables_bounds[i] = {lb, ub};
    }
    /* bounds of additional variables s, p and n are left unchanged */
    for (unsigned int i = problem.number_variables; i < current_iterate.x.size(); i++) {
        double lb = this->subproblem_variables_bounds[i].lb - current_iterate.x[i];
        double ub = this->subproblem_variables_bounds[i].ub - current_iterate.x[i];
        variables_bounds[i] = {lb, ub};
    }
    return variables_bounds;
}

void Sl1QP::recover_active_set(SubproblemSolution& solution, std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds) {
    ActiveSet active_set;
    // variables
    for (unsigned int i = 0; i < solution.x.size(); i++) {
        if (solution.x[i] == variables_bounds[i].lb) {
            active_set.at_lower_bound.push_back(i);
        }
        else if (solution.x[i] == variables_bounds[i].ub) {
            active_set.at_upper_bound.push_back(i);
        }
    }
    // constraints: only when p-n = 0
    for (unsigned int j = 0; j < solution.multipliers.constraints.size(); j++) {
        if (solution.x[this->positive_part_variables[j]] + solution.x[this->negative_part_variables[j]] == 0.) {
            active_set.at_lower_bound.push_back(solution.x.size() + j);
        }
    }
    solution.active_set = active_set;
    return;
}