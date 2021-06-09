#ifndef FEASIBILITYRESTORATION_H
#define FEASIBILITYRESTORATION_H

#include "FeasibilityStrategy.hpp"

class FeasibilityRestoration: public FeasibilityStrategy {
public:
   explicit FeasibilityRestoration(Subproblem& subproblem);
   std::vector<Direction> compute_feasible_directions(Problem& problem, Iterate& current_iterate, double trust_region_radius) override;
   double compute_predicted_reduction(Problem& problem, Iterate& current_iterate, Direction& direction, double step_length) override;
};

#endif //FEASIBILITYRESTORATION_H