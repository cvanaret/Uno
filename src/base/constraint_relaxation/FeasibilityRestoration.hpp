#ifndef FEASIBILITYRESTORATION_H
#define FEASIBILITYRESTORATION_H

#include "ConstraintRelaxationStrategy.hpp"

class FeasibilityRestoration : public ConstraintRelaxationStrategy {
public:
   explicit FeasibilityRestoration(Subproblem& subproblem);
   std::vector<Direction> compute_feasible_directions(Problem& problem, Iterate& current_iterate, double trust_region_radius) override;
   std::optional<Iterate> check_acceptance(Statistics& statistics, Problem& problem, Iterate& current_iterate, Direction& direction, double
   step_length) override;
   double compute_predicted_reduction(Problem& problem, Iterate& current_iterate, Direction& direction, double step_length) override;

private:
   Phase current_phase;
   //GlobalizationStrategy& phase_1_strategy;
   //GlobalizationStrategy& phase_2_strategy;
};

#endif //FEASIBILITYRESTORATION_H