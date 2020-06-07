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

Direction SLP::compute_step(Problem& problem, Iterate& current_iterate, double trust_region_radius) {
    Direction direction = this->compute_lp_step_(problem, this->solver, current_iterate, trust_region_radius);
    if (direction.status == INFEASIBLE) {
        /* infeasible subproblem during optimality phase */
        direction = this->restore_feasibility(problem, current_iterate, direction, trust_region_radius);
    }
    return direction;
}

Direction SLP::restore_feasibility(Problem& problem, Iterate& current_iterate, Direction& phase_2_direction, double trust_region_radius) {
   return this->compute_feasibility_lp_step_(problem, this->solver, current_iterate, phase_2_direction, trust_region_radius); 
}

/* private methods */

void SLP::evaluate_optimality_iterate_(Problem& problem, Iterate& current_iterate) {
    /* compute first-order information */
    current_iterate.compute_objective_gradient(problem);
    current_iterate.compute_constraints_jacobian(problem);
}

void SLP::evaluate_feasibility_iterate_(Problem& problem, Iterate& current_iterate, Direction& /*phase_2_direction*/) {
    /* compute first-order information */
    current_iterate.compute_constraints_jacobian(problem);
}