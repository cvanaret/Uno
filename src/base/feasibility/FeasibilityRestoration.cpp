#include "FeasibilityRestoration.hpp"

FeasibilityRestoration::FeasibilityRestoration(Subproblem& subproblem) : FeasibilityStrategy(subproblem) {
}

std::vector<Direction> FeasibilityRestoration::compute_feasible_directions(Problem& /*problem*/, Iterate& /*current_iterate*/, double
/*trust_region_radius*/) {
   return std::vector<Direction>();
}

double FeasibilityRestoration::compute_predicted_reduction(Problem& /*problem*/, Iterate& /*current_iterate*/, Direction& /*direction*/,
      double /*step_length*/) {
   return 0.;
}