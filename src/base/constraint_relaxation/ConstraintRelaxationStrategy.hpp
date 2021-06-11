#ifndef CONSTRAINTRELAXATIONSTRATEGY_H
#define CONSTRAINTRELAXATIONSTRATEGY_H

#include <vector>
#include <cmath>
#include "Statistics.hpp"
#include "Subproblem.hpp"
#include "Direction.hpp"
#include "Problem.hpp"
#include "Iterate.hpp"

class ConstraintRelaxationStrategy {
public:
   Subproblem& subproblem;

   explicit ConstraintRelaxationStrategy(Subproblem& subproblem);
   virtual Iterate initialize(Statistics& statistics, const Problem& problem, std::vector<double>& x, Multipliers& multipliers) = 0;
   virtual Direction compute_feasible_direction(const Problem& problem, Iterate& current_iterate, double trust_region_radius) = 0;
   virtual std::optional<Iterate> check_acceptance(Statistics& statistics, const Problem& problem, Iterate& current_iterate, Direction& direction,
         double step_length) = 0;
   virtual double compute_predicted_reduction(const Problem& problem, Iterate& current_iterate, Direction& direction, double step_length) = 0;

protected:
   // preallocated vector to receive the trial primal variables
   std::vector<double> trial_primals_;
};

#endif //CONSTRAINTRELAXATIONSTRATEGY_H
