#ifndef UNO_L1RELAXATION_H
#define UNO_L1RELAXATION_H

#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/strategy/GlobalizationStrategy.hpp"
#include "tools/Options.hpp"

struct l1RelaxationParameters {
   bool fixed_parameter;
   double decrease_factor;
   double epsilon1;
   double epsilon2;
   double small_threshold;
};

class l1Relaxation : public ConstraintRelaxationStrategy {
public:
   l1Relaxation(const Problem& problem, const Options& options);
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

protected:
   const std::unique_ptr<GlobalizationStrategy> globalization_strategy;
   double penalty_parameter;
   const l1RelaxationParameters parameters;
   // preallocated temporary multipliers
   std::vector<double> constraint_multipliers;

   static double compute_predicted_reduction(const Problem& problem, Iterate& current_iterate, const Direction& direction,
         PredictedReductionModel& predicted_reduction_model, double step_length);
   Direction solve_subproblem(Statistics& statistics, Iterate& current_iterate, double current_penalty_parameter);
   Direction resolve_subproblem(Statistics& statistics, Iterate& current_iterate, double current_penalty_parameter);
   Direction solve_with_steering_rule(Statistics& statistics, Iterate& current_iterate);
   void decrease_parameter_aggressively(Iterate& current_iterate, const Direction& direction_lowest_violation);
   [[nodiscard]] bool linearized_residual_sufficient_decrease(const Iterate& current_iterate, double linearized_residual, double residual_lowest_violation) const;
   [[nodiscard]] bool objective_sufficient_decrease(const Iterate& current_iterate, const Direction& direction, const Direction& direction_lowest_violation) const;
   double compute_error(Iterate& current_iterate, const Multipliers& multiplier_displacements, double current_penalty_parameter);
   void set_multipliers(const Iterate& current_iterate, std::vector<double>& constraint_multipliers);
};

#endif //UNO_L1RELAXATION_H