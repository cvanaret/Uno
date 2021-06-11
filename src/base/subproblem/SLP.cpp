#include <map>
#include "SLP.hpp"
#include "QPSolverFactory.hpp"

SLP::SLP(const Problem& problem, std::string QP_solver_name, bool /*use_trust_region*/, bool scale_residuals):
ActiveSetMethod(problem, scale_residuals),
solver(QPSolverFactory::create(QP_solver_name, problem.number_variables, problem.number_constraints, 0, false)) {
}

void SLP::generate(const Problem& /*problem*/, const Iterate& /*current_iterate*/, double /*objective_multiplier*/, double
/*trust_region_radius*/) {
}

void SLP::update_objective_multipliers(const Problem& /*problem*/, const Iterate& /*current_iterate*/, double /*objective_multiplier*/) {
}

Direction SLP::compute_direction(const Problem& problem, Iterate& current_iterate, double trust_region_radius) {
    Direction direction = this->compute_lp_step_(problem, *this->solver, current_iterate, trust_region_radius);
    if (direction.status != INFEASIBLE) {
       return direction;
    }
    else {
        /* infeasible subproblem during optimality phase */
        return this->restore_feasibility(problem, current_iterate, direction, trust_region_radius);
    }
}

Direction SLP::restore_feasibility(const Problem& problem, Iterate& current_iterate, Direction& phase_2_direction, double trust_region_radius) {
   Direction direction = this->compute_l1lp_step_(problem, *this->solver, current_iterate, phase_2_direction, trust_region_radius);
   direction.is_relaxed = true;
   return direction;
}

/* private methods */

void SLP::evaluate_optimality_iterate_(const Problem& problem, Iterate& current_iterate) {
    /* compute first-order information */
    current_iterate.compute_objective_gradient(problem);
    current_iterate.compute_constraints_jacobian(problem);
}

void SLP::evaluate_feasibility_iterate_(const Problem& problem, Iterate& current_iterate, Direction& /*phase_2_direction*/) {
    /* compute first-order information */
    current_iterate.compute_constraints_jacobian(problem);
}

