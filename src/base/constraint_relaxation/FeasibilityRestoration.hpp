#ifndef FEASIBILITYRESTORATION_H
#define FEASIBILITYRESTORATION_H

#include "ConstraintRelaxationStrategy.hpp"
#include "GlobalizationStrategy.hpp"

enum Phase {FEASIBILITY_RESTORATION = 1, OPTIMALITY = 2};

class FeasibilityRestoration : public ConstraintRelaxationStrategy {
public:
   explicit FeasibilityRestoration(const std::string& constraint_relaxation_strategy, Subproblem& subproblem, const std::map<std::string, std::string>& options);
   Iterate initialize(Statistics& statistics, const Problem& problem, std::vector<double>& x, Multipliers& multipliers) override;
   Direction compute_feasible_direction(const Problem& problem, Iterate& current_iterate, double trust_region_radius) override;
   bool is_acceptable(Statistics& statistics, const Problem& problem, Iterate& current_iterate, Iterate& trial_iterate, Direction& direction, double
   step_length) override;
   double compute_predicted_reduction(const Problem& problem, Iterate& current_iterate, Direction& direction, double step_length) override;

private:
   std::unique_ptr<GlobalizationStrategy> phase_1_strategy;
   std::unique_ptr<GlobalizationStrategy> phase_2_strategy;
   Phase current_phase;

   static void update_restoration_multipliers_(Iterate& trial_iterate, const ConstraintPartition& constraint_partition);
};

#endif //FEASIBILITYRESTORATION_H