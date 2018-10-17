#include <cmath>
#include "TwoPhaseStrategy.hpp"
#include "Logger.hpp"

TwoPhaseStrategy::TwoPhaseStrategy(Subproblem& subproblem, TwoPhaseConstants& constants, double tolerance) :
GlobalizationStrategy(subproblem, tolerance), phase(OPTIMALITY), constants(constants) {
}

LocalSolution TwoPhaseStrategy::compute_step(Problem& problem, Iterate& current_iterate, double radius) {
    //double objective_multiplier = (phase == OPTIMALITY) ? problem.obj_sign : 0.;
    LocalSolution solution = this->subproblem.compute_optimality_step(problem, current_iterate, radius);
    solution.phase = OPTIMALITY;
    DEBUG << solution;

    if (solution.phase_1_required) {
        /* infeasible subproblem during optimality phase */
        DEBUG << "Moving on to the feasibility step\n";

        /* different Hessian, but same constraint Jacobian */
        current_iterate.is_hessian_computed = false;
        /* partition of feasible and infeasible phase II constraints */
        ConstraintPartition constraint_partition = solution.constraint_partition;

        /* compute the step in phase 1, starting from infeasible solution */
        solution = this->subproblem.compute_infeasibility_step(problem, current_iterate, radius, solution);
        solution.constraint_partition = constraint_partition;
        solution.phase = RESTORATION;
        DEBUG << solution;
    }
    /* from this point on, the step is feasible */
    return solution;
}

void TwoPhaseStrategy::update_restoration_multipliers(Iterate& trial_iterate, ConstraintPartition& constraint_partition) {
    for (unsigned int k = 0; k < constraint_partition.infeasible_set.size(); k++) {
        int j = constraint_partition.infeasible_set[k];
        if (constraint_partition.constraint_status[j] == INFEASIBLE_UPPER) {
            trial_iterate.constraint_multipliers[j] = -1.;
        }
        else {
            trial_iterate.constraint_multipliers[j] = 1.;
        }
    }
    return;
}