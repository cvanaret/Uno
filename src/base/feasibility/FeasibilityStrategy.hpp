#ifndef INFEASIBILITYSTRATEGY_H
#define INFEASIBILITYSTRATEGY_H

#include <vector>
#include <cmath>
#include "Subproblem.hpp"
#include "Direction.hpp"
#include "Problem.hpp"
#include "Iterate.hpp"

class FeasibilityStrategy {
public:
   Subproblem& subproblem;

   explicit FeasibilityStrategy(Subproblem& subproblem);
   virtual std::vector<Direction> compute_feasible_directions(Problem& problem, Iterate& current_iterate, double trust_region_radius) = 0;
   virtual double compute_predicted_reduction(Problem& problem, Iterate& current_iterate, Direction& direction, double step_length) = 0;
};

#endif //INFEASIBILITYSTRATEGY_H
