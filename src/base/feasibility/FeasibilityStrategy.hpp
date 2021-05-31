#ifndef INFEASIBILITYSTRATEGY_H
#define INFEASIBILITYSTRATEGY_H

#include <vector>
#include <cmath>
#include "Subproblem.hpp"
#include "Direction.hpp"
#include "Problem.hpp"
#include "Iterate.hpp"

class InfeasibilityStrategy {
public:
   Subproblem& subproblem;

   explicit InfeasibilityStrategy(Subproblem& subproblem);
   virtual std::vector<Direction> compute_directions(Problem& problem, Iterate& current_iterate, double trust_region_radius) = 0;
};

#endif //INFEASIBILITYSTRATEGY_H
