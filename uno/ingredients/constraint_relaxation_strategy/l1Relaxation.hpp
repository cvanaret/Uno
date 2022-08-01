// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_L1RELAXATION_H
#define UNO_L1RELAXATION_H

#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/globalization_strategy/GlobalizationStrategy.hpp"
#include "reformulation/OptimalityProblem.hpp"
#include "reformulation/l1RelaxedProblem.hpp"
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
   l1Relaxation(const Model& model, const Options& options);
   void initialize(Statistics& statistics, Iterate& first_iterate) override;

   void set_variable_bounds(const Iterate& current_iterate, double trust_region_radius) override;

   // direction computation
   [[nodiscard]] Direction compute_feasible_direction(Statistics& statistics, Iterate& current_iterate) override;
   [[nodiscard]] Direction solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate) override;
   [[nodiscard]] Direction solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate,
         const std::vector<double>& initial_point) override;
   [[nodiscard]] Direction compute_second_order_correction(Iterate& trial_iterate) override;

   // trial iterate acceptance
   [[nodiscard]] bool is_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         PredictedOptimalityReductionModel& predicted_optimality_reduction_model, double step_length) override;
   [[nodiscard]] PredictedOptimalityReductionModel generate_predicted_optimality_reduction_model(const Direction& direction) const override;
   void register_accepted_iterate(Iterate& iterate) override;

   [[nodiscard]] size_t get_hessian_evaluation_count() const override;
   [[nodiscard]] size_t get_number_subproblems_solved() const override;

protected:
   const OptimalityProblem optimality_problem;
   l1RelaxedProblem relaxed_problem;
   std::unique_ptr<Subproblem> subproblem;
   const std::unique_ptr<GlobalizationStrategy> globalization_strategy;
   double penalty_parameter;
   const l1RelaxationParameters parameters;
   // preallocated temporary multipliers
   std::vector<double> constraint_multipliers;
   std::vector<double> lower_bound_multipliers;
   std::vector<double> upper_bound_multipliers;
   // statistics table
   int statistics_penalty_parameter_column_order;

   Direction solve_subproblem(Statistics& statistics, Iterate& current_iterate, double current_penalty_parameter);
   Direction solve_with_steering_rule(Statistics& statistics, Iterate& current_iterate);
   void decrease_parameter_aggressively(Iterate& current_iterate, const Direction& direction_lowest_violation);
   [[nodiscard]] bool linearized_residual_sufficient_decrease(const Iterate& current_iterate, double linearized_residual, double residual_lowest_violation) const;
   [[nodiscard]] bool objective_sufficient_decrease(const Iterate& current_iterate, const Direction& direction, const Direction& direction_lowest_violation) const;
   double compute_error(Iterate& current_iterate, const Multipliers& multiplier_displacements);
   void set_multipliers(const Iterate& current_iterate, std::vector<double>& constraint_multipliers);
   double compute_infeasibility_measure(Iterate& iterate) override;
};

#endif //UNO_L1RELAXATION_H