#include "FeasibilityRestoration.hpp"

FeasibilityRestoration::FeasibilityRestoration(Subproblem& subproblem) : FeasibilityStrategy(subproblem) {
}

std::vector<Direction> FeasibilityRestoration::compute_feasible_directions(Problem& problem, Iterate& current_iterate, double
trust_region_radius) {
   std::vector<Direction> directions = this->subproblem.compute_directions(problem, current_iterate, trust_region_radius);
   // TODO handle multiple directions
   Direction& direction = directions[0];

   if (direction.status != INFEASIBLE) {
      direction.phase = OPTIMALITY;
      direction.objective_multiplier = problem.objective_sign;
      return directions;
   }
   else {
      /* infeasible subproblem: switch to restoration phase */
      return this->subproblem.restore_feasibility(problem, current_iterate, direction, trust_region_radius);
   }
}

double FeasibilityRestoration::compute_predicted_reduction(Problem& /*problem*/, Iterate& /*current_iterate*/, Direction& direction,
      double step_length) {
   return direction.predicted_reduction(step_length);
}