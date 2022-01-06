#ifndef UNO_FEASIBILITYRESTORATION_H
#define UNO_FEASIBILITYRESTORATION_H

#include <optional>
#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/strategy/GlobalizationStrategy.hpp"
#include "tools/Options.hpp"

enum Phase {FEASIBILITY_RESTORATION = 1, OPTIMALITY = 2};

class FeasibilityRestoration : public ConstraintRelaxationStrategy {
public:
   FeasibilityRestoration(const Problem& problem, const Options& options);
   void initialize(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& first_iterate) override;

   // direction computation
   void create_current_subproblem(const Problem& problem, const Scaling& scaling, Iterate& current_iterate, double trust_region_radius) override;
   Direction compute_feasible_direction(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& current_iterate) override;
   Direction solve_feasibility_problem(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& current_iterate,
         const Direction& direction) override;

   bool is_acceptable(Statistics& statistics, const Problem& problem, const Scaling& scaling, Iterate& current_iterate, Iterate& trial_iterate,
         const Direction& direction, PredictedReductionModel& predicted_reduction_model, double step_length) override;

private:
   const std::unique_ptr<GlobalizationStrategy> phase_1_strategy;
   const std::unique_ptr<GlobalizationStrategy> phase_2_strategy;
   Phase current_phase{OPTIMALITY};

   void form_feasibility_problem(const Problem& problem, const Scaling& scaling, Iterate& current_iterate,
         const std::vector<double>& phase_2_primal_direction, const std::optional<ConstraintPartition>& optional_constraint_partition);
   GlobalizationStrategy& switch_phase(const Problem& problem, const Scaling& scaling, Iterate& current_iterate, const Direction& direction);
   static void set_restoration_multipliers(std::vector<double>& constraints_multipliers, const ConstraintPartition& constraint_partition);
   void compute_infeasibility_measures(const Problem& problem, const Scaling& scaling, Iterate& iterate,
         const std::optional<ConstraintPartition>& optional_constraint_partition);
};

#endif //UNO_FEASIBILITYRESTORATION_H