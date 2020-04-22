#include <cmath>
#include <map>
#include "SLP.hpp"
#include "Constraint.hpp"
#include "Utils.hpp"
#include "Logger.hpp"

SLP::SLP(QPSolver& solver) : ActiveSetMethod(solver) {
}

double SLP::compute_predicted_reduction(Problem& problem, Iterate& current_iterate, SubproblemSolution& solution, double step_length) {
    if (step_length == 1.) {
        return -solution.objective;
    }
    else {
        double linear_term = dot(solution.x, current_iterate.objective_gradient);
        return -step_length*linear_term;
    }
}

bool SLP::phase_1_required(SubproblemSolution& solution) {
    return (solution.status == INFEASIBLE);
}

/* private methods */

void SLP::evaluate_optimality_iterate(Problem& problem, Iterate& current_iterate) {
    /* compute first-order information */
    current_iterate.compute_objective_gradient(problem);
    current_iterate.compute_constraints_jacobian(problem);
}

void SLP::evaluate_feasibility_iterate(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution) {
    /* compute first-order information */
    current_iterate.compute_constraints_jacobian(problem);
}

SubproblemSolution SLP::solve_subproblem(std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds, Iterate& current_iterate, std::vector<double>& d0) {
    return this->solver.solve_LP(variables_bounds, constraints_bounds, current_iterate.objective_gradient, current_iterate.constraints_jacobian, d0);
}