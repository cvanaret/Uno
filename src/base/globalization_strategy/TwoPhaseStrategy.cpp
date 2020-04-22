#include <cmath>
#include "TwoPhaseStrategy.hpp"
#include "Logger.hpp"

TwoPhaseStrategy::TwoPhaseStrategy(Subproblem& subproblem, TwoPhaseConstants& constants, double tolerance) :
GlobalizationStrategy(subproblem, tolerance), phase(OPTIMALITY), constants(constants) {
}

SubproblemSolution TwoPhaseStrategy::compute_step(Problem& problem, Iterate& current_iterate, double trust_region_radius) {
    //double objective_multiplier = (phase == OPTIMALITY) ? problem.obj_sign : 0.;
    SubproblemSolution solution = this->subproblem.compute_optimality_step(problem, current_iterate, trust_region_radius);
    solution.phase = OPTIMALITY;
    DEBUG << solution;

    if (solution.phase_1_required) {
        /* infeasible subproblem during optimality phase */
        solution = this->restore_feasibility(problem, current_iterate, solution, trust_region_radius);
    }
    /* from this point on, the step is feasible */
    return solution;
}

SubproblemSolution TwoPhaseStrategy::restore_feasibility(Problem& problem, Iterate& current_iterate, SubproblemSolution& phase_II_solution, double trust_region_radius) {
    /* different Hessian, but same constraint Jacobian */
    current_iterate.is_hessian_computed = false;
    /* partition of feasible and infeasible phase II constraints */
    //ConstraintPartition constraint_partition = phase_II_solution.constraint_partition;

    /* compute the step in phase 1, starting from infeasible solution */
    DEBUG << "Creating the restoration problem\n";
    SubproblemSolution solution = this->subproblem.compute_infeasibility_step(problem, current_iterate, phase_II_solution, trust_region_radius);
    solution.phase = RESTORATION;
    DEBUG << solution;
    return solution;
}

void TwoPhaseStrategy::update_restoration_multipliers(Iterate& trial_iterate, ConstraintPartition& constraint_partition) {
    for (int j: constraint_partition.infeasible) {
        if (constraint_partition.constraint_feasibility[j] == INFEASIBLE_UPPER) {
            trial_iterate.multipliers.constraints[j] = -1.;
        }
        else {
            trial_iterate.multipliers.constraints[j] = 1.;
        }
    }
    return;
}
