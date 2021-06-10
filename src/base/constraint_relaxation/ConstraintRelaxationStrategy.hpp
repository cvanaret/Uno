#ifndef CONSTRAINTRELAXATIONSTRATEGY_H
#define CONSTRAINTRELAXATIONSTRATEGY_H

#include <vector>
#include <cmath>
#include "Subproblem.hpp"
#include "Direction.hpp"
#include "Problem.hpp"
#include "Iterate.hpp"

class ConstraintRelaxationStrategy {
public:
   Subproblem& subproblem;

   explicit ConstraintRelaxationStrategy(Subproblem& subproblem);
   virtual std::vector<Direction> compute_feasible_directions(Problem& problem, Iterate& current_iterate, double trust_region_radius) = 0;
   virtual double compute_predicted_reduction(Problem& problem, Iterate& current_iterate, Direction& direction, double step_length) = 0;
};

#endif //CONSTRAINTRELAXATIONSTRATEGY_H
