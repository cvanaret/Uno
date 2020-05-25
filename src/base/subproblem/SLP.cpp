#include <cmath>
#include <map>
#include "SLP.hpp"
#include "Constraint.hpp"
#include "Utils.hpp"
#include "Logger.hpp"
#include "QPSolverFactory.hpp"

SLP::SLP(Problem& problem, std::string QP_solver_name, bool /*use_trust_region*/, bool scale_residuals):
ActiveSetMethod(problem, scale_residuals),
solver(QPSolverFactory::create(QP_solver_name, problem.number_variables, problem.number_constraints, 0, false)) {
}

SubproblemSolution SLP::compute_step(Problem& problem, Iterate& current_iterate, double trust_region_radius) {
    SubproblemSolution solution = this->compute_lp_step_(problem, this->solver, current_iterate, trust_region_radius);
    if (solution.status == INFEASIBLE) {
        /* infeasible subproblem during optimality phase */
        solution = this->restore_feasibility(problem, current_iterate, solution, trust_region_radius);
    }
    return solution;
}

SubproblemSolution SLP::restore_feasibility(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution, double trust_region_radius) {
   return this->compute_feasibility_lp_step_(problem, this->solver, current_iterate, phase_II_solution, trust_region_radius); 
}

/* private methods */

void SLP::evaluate_optimality_iterate_(Problem& problem, Iterate& current_iterate) {
    /* compute first-order information */
    current_iterate.compute_objective_gradient(problem);
    current_iterate.compute_constraints_jacobian(problem);
}

void SLP::evaluate_feasibility_iterate_(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution) {
    /* compute first-order information */
    current_iterate.compute_constraints_jacobian(problem);
}