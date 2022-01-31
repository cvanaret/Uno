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
   void initialize(Statistics& statistics, Iterate& first_iterate) override;

   // direction computation
   void create_current_subproblem(Iterate& current_iterate, double trust_region_radius) override;
   [[nodiscard]] Direction compute_feasible_direction(Statistics& statistics, Iterate& current_iterate) override;
   [[nodiscard]] Direction solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate,
         const std::optional<std::vector<double>>& optional_phase_2_solution) override;
   [[nodiscard]] Direction compute_second_order_correction(Iterate& trial_iterate) override;

   [[nodiscard]] bool is_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         PredictedReductionModel& predicted_reduction_model, double step_length) override;
   [[nodiscard]] PredictedReductionModel generate_predicted_reduction_model(const Iterate& current_iterate, const Direction& direction) const override;
   void register_accepted_iterate(Iterate& iterate) override;

private:
   const std::unique_ptr<GlobalizationStrategy> phase_1_strategy;
   const std::unique_ptr<GlobalizationStrategy> phase_2_strategy;
   Phase current_phase{OPTIMALITY};

   void create_current_feasibility_problem(Iterate& current_iterate, const std::optional<std::vector<double>>& optional_phase_2_solution);
   GlobalizationStrategy& switch_phase(Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction);
};

#endif //UNO_FEASIBILITYRESTORATION_H