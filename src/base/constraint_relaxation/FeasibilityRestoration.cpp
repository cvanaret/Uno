#include "FeasibilityRestoration.hpp"

FeasibilityRestoration::FeasibilityRestoration(Subproblem& subproblem) : ConstraintRelaxationStrategy(subproblem), current_phase(OPTIMALITY) {
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
      directions = this->subproblem.restore_feasibility(problem, current_iterate, direction, trust_region_radius);
      /* infeasible subproblem: go from phase II (optimality) to I (restoration) */
      DEBUG << "Switching from optimality to restoration phase\n";
      this->current_phase = FEASIBILITY_RESTORATION;
      return directions;
   }
}

std::optional<Iterate> FeasibilityRestoration::check_acceptance(Statistics& statistics, Problem& problem, Iterate& current_iterate, Direction& direction,
      double step_length) {
   return std::optional<Iterate>();
}

double FeasibilityRestoration::compute_predicted_reduction(Problem& /*problem*/, Iterate& /*current_iterate*/, Direction& direction,
      double step_length) {
   return direction.predicted_reduction(step_length);
}