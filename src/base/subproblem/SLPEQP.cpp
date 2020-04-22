#include <cmath>
#include <map>
#include "SLPEQP.hpp"
#include "SLP.hpp"
#include "SQP.hpp"
#include "Constraint.hpp"
#include "Utils.hpp"
#include "Logger.hpp"

SLPEQP::SLPEQP(QPSolver& solver, HessianEvaluation& hessian_evaluation) :
Subproblem(1), solver(solver), lp_subproblem(SLP(solver)), eqp_subproblem(SQP(solver, hessian_evaluation)) {
}

Iterate SLPEQP::initialize(Problem& problem, std::vector<double>& x, Multipliers& multipliers, int number_variables, bool use_trust_region) {
    // register the original bounds
    this->subproblem_variables_bounds = problem.variables_bounds;
    
    Iterate first_iterate(problem, x, multipliers, this->residual_norm);
    /* compute the optimality and feasibility measures of the initial point */
    this->compute_measures(problem, first_iterate);

    /* allocate the QP solver */
    this->solver.allocate(number_variables, problem.number_constraints);
    return first_iterate;
}

SubproblemSolution SLPEQP::compute_optimality_step(Problem& problem, Iterate& current_iterate, double trust_region_radius) {
    /* compute first-order information */
    current_iterate.compute_objective_gradient(problem);
    current_iterate.compute_constraints_jacobian(problem);
    
    /***********/
    /* LP part */
    /***********/
    
    //std::cout << "SOLVING SLPEQP.LP\n";
    
    /* bounds of the variables */
    std::vector<Range> lp_variables_bounds = this->generate_variables_bounds(current_iterate, trust_region_radius);
    
    /* bounds of the linearized constraints */
    std::vector<Range> constraints_bounds = Subproblem::generate_constraints_bounds(problem, current_iterate.constraints);

    /* generate the initial point */
    std::vector<double> d0(lp_variables_bounds.size()); // = {0.}

    /* solve the LP */
    SubproblemSolution solution_LP = this->solver.solve_LP(lp_variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate.constraints_jacobian, d0);
    solution_LP.phase_1_required = this->phase_1_required(solution_LP);
    print_vector(std::cout, solution_LP.x);
    
    /************/
    /* EQP part */
    /************/
    
    //std::cout << "SOLVING SLPEQP.EQP\n";
    /* compute second-order information */
    current_iterate.compute_hessian(problem, problem.objective_sign, current_iterate.multipliers.constraints);
    
    /* fix active constraints */
    // TODO remove inactive constraints!!
    // TODO reduced Hessian
    this->fix_active_constraints(problem, solution_LP.active_set, lp_variables_bounds, constraints_bounds);
    
    /* generate the initial point */
    d0 = solution_LP.x;
    
    /* solve the EQP */
    SubproblemSolution solution_EQP = this->solver.solve_QP(problem.variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate.constraints_jacobian, current_iterate.hessian, d0);
    solution_EQP.phase_1_required = this->phase_1_required(solution_EQP);
    //print_vector(std::cout, solution_EQP.x);
    this->number_subproblems_solved++;
    
    return solution_EQP;
}

SubproblemSolution SLPEQP::compute_infeasibility_step(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution, double trust_region_radius) {
    // TODO
    return this->compute_optimality_step(problem, current_iterate, trust_region_radius);
}

double SLPEQP::compute_predicted_reduction(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution, double step_length) {
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

void SLPEQP::compute_measures(Problem& problem, Iterate& iterate) {
    iterate.feasibility_measure = iterate.residual;
    iterate.optimality_measure = iterate.objective;
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