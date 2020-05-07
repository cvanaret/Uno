#include <cmath>
#include <map>
#include "SLPEQP.hpp"
#include "SLP.hpp"
#include "SQP.hpp"
#include "Constraint.hpp"
#include "Utils.hpp"
#include "Logger.hpp"

SLPEQP::SLPEQP(Problem& problem, std::string QP_solver_name, std::string hessian_evaluation_method):
Subproblem("l1"),
lp_subproblem(SLP(problem, QP_solver_name)),
eqp_subproblem(SQP(problem, QP_solver_name, hessian_evaluation_method)) {
}

Iterate SLPEQP::initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers, bool use_trust_region) {
    // register the original bounds
    this->subproblem_variables_bounds = problem.variables_bounds;
    
    this->lp_subproblem.initialize(problem, x, multipliers, use_trust_region);
    this->eqp_subproblem.initialize(problem, x, multipliers, use_trust_region);

    Iterate first_iterate(x, multipliers);
    /* compute the optimality and feasibility measures of the initial point */
    this->compute_optimality_measures(problem, first_iterate);
    
    return first_iterate;
}

SubproblemSolution SLPEQP::compute_optimality_step(Problem& problem, Iterate& current_iterate, double trust_region_radius) {
    /* compute first-order information */
    current_iterate.compute_objective_gradient(problem);
    current_iterate.compute_constraints(problem);
    current_iterate.compute_constraints_jacobian(problem);

    /***********/
    /* LP part */
    /***********/
    DEBUG << "SOLVING SLPEQP.LP\n";

    /* solve the LP */
    SubproblemSolution solution_LP = this->lp_subproblem.compute_optimality_step(problem, current_iterate, trust_region_radius);
    solution_LP.phase_1_required = this->phase_1_required(solution_LP);
    print_vector(DEBUG, solution_LP.x);

    /************/
    /* EQP part */
    /************/
    DEBUG << "SOLVING SLPEQP.EQP\n";

    /* fix active constraints */
    // TODO remove inactive constraints!!
    // TODO reduced Hessian
    this->fix_active_constraints(problem, solution_LP.active_set, lp_variables_bounds, constraints_bounds);

    /* generate the initial point */
    d0 = solution_LP.x;

    /* solve the EQP */
    SubproblemSolution solution_EQP = this->solver.solve_QP(problem.variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate.constraints_jacobian, current_iterate.hessian, d0);
    solution_EQP.phase_1_required = this->phase_1_required(solution_EQP);
    //print_vector(DEBUG, solution_EQP.x);
    this->number_subproblems_solved++;

    return solution_EQP;
}

SubproblemSolution SLPEQP::compute_infeasibility_step(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution, double trust_region_radius) {
    // TODO
    return this->compute_optimality_step(problem, current_iterate, trust_region_radius);
}

double SLPEQP::compute_predicted_reduction(Problem& /*problem*/, Iterate& /*current_iterate*/, SubproblemSolution& solution, double step_length) {
    return -solution.objective;
}

void SLPEQP::compute_optimality_measures(Problem& problem, Iterate& iterate) {
    // feasibility
    iterate.compute_constraint_residual(problem, this->residual_norm);
    iterate.feasibility_measure = iterate.constraint_residual;
    // optimality
    iterate.compute_objective(problem);
    iterate.optimality_measure = iterate.objective;
    return;
}

void SLPEQP::compute_infeasibility_measures(Problem& problem, Iterate& iterate, SubproblemSolution& solution) {
    iterate.compute_constraints(problem);
    iterate.feasibility_measure = problem.compute_constraint_residual(solution.constraint_partition, iterate.constraints, this->residual_norm);
    iterate.optimality_measure = problem.compute_constraint_residual(solution.constraint_partition, iterate.constraints, this->residual_norm);
    return;
}

bool SLPEQP::phase_1_required(SubproblemSolution& solution) {
    return solution.status == INFEASIBLE;
}

/* private methods */

void SLPEQP::fix_active_constraints(Problem& problem, ActiveSet& active_set, std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds) {
    /* active lower bounds */
    for (int i: active_set.at_lower_bound) {
        if (i < problem.number_variables) {
            if (problem.variables_bounds[i].lb == variables_bounds[i].lb) { // bound constraint
                variables_bounds[i].ub = variables_bounds[i].lb;
            }
            // otherwise, trust region
        }
        else { // general constraint
            int j = i - problem.number_variables;
            constraints_bounds[j].ub = constraints_bounds[j].lb;
        }
    }
    /* active upper bounds */
    for (int i: active_set.at_upper_bound) {
        if (i < problem.number_variables) {
            if (problem.variables_bounds[i].ub == variables_bounds[i].ub) { // bound constraint
                variables_bounds[i].lb = variables_bounds[i].ub;
            }
            // otherwise, trust region
        }
        else { // general constraint
            int j = i - problem.number_variables;
            constraints_bounds[j].lb = constraints_bounds[j].ub;
        }
    }
    return;
}