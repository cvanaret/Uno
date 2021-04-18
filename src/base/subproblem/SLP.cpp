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

std::vector<Direction> SLP::compute_directions(Problem& problem, Iterate& current_iterate, double trust_region_radius) {
    Direction direction = this->compute_lp_step_(problem, *this->solver, current_iterate, trust_region_radius);
    if (direction.status == INFEASIBLE) {
        /* infeasible subproblem during optimality phase */
        return this->restore_feasibility(problem, current_iterate, direction, trust_region_radius);
    }
    else {
        direction.phase = OPTIMALITY;
        return std::vector<Direction>{direction};
    }
}

std::vector<Direction> SLP::restore_feasibility(Problem& problem, Iterate& current_iterate, Direction& phase_2_direction, double trust_region_radius) {
   Direction direction = this->compute_l1lp_step_(problem, *this->solver, current_iterate, phase_2_direction, trust_region_radius);
   direction.phase = RESTORATION;
   return std::vector<Direction>{direction};
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
