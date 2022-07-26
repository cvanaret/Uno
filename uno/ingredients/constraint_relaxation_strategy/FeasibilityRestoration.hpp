// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FEASIBILITYRESTORATION_H
#define UNO_FEASIBILITYRESTORATION_H

#include <optional>
#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/globalization_strategy/GlobalizationStrategy.hpp"
#include "ingredients/constraint_relaxation_strategy/OptimalityProblem.hpp"
#include "ingredients/constraint_relaxation_strategy/l1RelaxedProblem.hpp"
#include "tools/Options.hpp"

enum Phase {FEASIBILITY_RESTORATION = 1, OPTIMALITY = 2};

class FeasibilityRestoration : public ConstraintRelaxationStrategy {
public:
   FeasibilityRestoration(const Model& model, const Options& options);
   void initialize(Statistics& statistics, Iterate& first_iterate) override;

   void set_variable_bounds(const Iterate& current_iterate, double trust_region_radius) override;

   // direction computation
   [[nodiscard]] Direction compute_feasible_direction(Statistics& statistics, Iterate& current_iterate) override;
   [[nodiscard]] Direction solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate,
         const std::optional<std::vector<double>>& optional_phase_2_solution) override;
   [[nodiscard]] Direction compute_second_order_correction(Iterate& trial_iterate) override;

   [[nodiscard]] bool is_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         PredictedOptimalityReductionModel& predicted_optimality_reduction_model, double step_length) override;
   [[nodiscard]] PredictedOptimalityReductionModel generate_predicted_optimality_reduction_model(const Direction& direction) const override;
   void register_accepted_iterate(Iterate& iterate) override;

   [[nodiscard]] size_t get_hessian_evaluation_count() const override;
   [[nodiscard]] size_t get_number_subproblems_solved() const override;

private:
   const OptimalityProblem optimality_problem;
   l1RelaxedProblem feasibility_problem;
   std::unique_ptr<Subproblem> subproblem;
   const std::unique_ptr<GlobalizationStrategy> phase_1_strategy;
   const std::unique_ptr<GlobalizationStrategy> phase_2_strategy;
   Phase current_phase{OPTIMALITY};
   // statistics table
   int statistics_restoration_phase_column_order;

   [[nodiscard]] const NonlinearProblem& get_current_reformulated_problem() const;
   [[nodiscard]] Direction solve_optimality_problem(Statistics& statistics, Iterate& current_iterate);
   [[nodiscard]] GlobalizationStrategy& switch_phase(Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction);
   void switch_to_feasibility_restoration(Iterate& current_iterate, const std::vector<size_t>& infeasible_constraints);
   void switch_to_optimality(Iterate& current_iterate, Iterate& trial_iterate);
   [[nodiscard]] double compute_infeasibility_measure(Iterate& iterate) override;
   [[nodiscard]] double compute_optimality_measure(Iterate& iterate, const std::vector<size_t>& infeasible_constraints);
};

#endif //UNO_FEASIBILITYRESTORATION_H