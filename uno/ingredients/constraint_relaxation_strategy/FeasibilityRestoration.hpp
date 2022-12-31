// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FEASIBILITYRESTORATION_H
#define UNO_FEASIBILITYRESTORATION_H

#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/globalization_strategy/GlobalizationStrategy.hpp"
#include "reformulation/OptimalityProblem.hpp"
#include "reformulation/l1RelaxedProblem.hpp"
#include "tools/Options.hpp"

enum class Phase {FEASIBILITY_RESTORATION = 1, OPTIMALITY = 2};

class FeasibilityRestoration : public ConstraintRelaxationStrategy {
public:
   FeasibilityRestoration(const Model& model, const Options& options);
   void initialize(Statistics& statistics, Iterate& first_iterate) override;

   void set_trust_region_radius(double trust_region_radius) override;

   // direction computation
   [[nodiscard]] Direction compute_feasible_direction(Statistics& statistics, Iterate& current_iterate) override;
   [[nodiscard]] Direction solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate) override;
   [[nodiscard]] Direction solve_feasibility_problem(Statistics& statistics, Iterate& current_iterate, const std::vector<double>& initial_point)
      override;
   [[nodiscard]] Direction compute_second_order_correction(Iterate& trial_iterate) override;

   // trial iterate acceptance
   void compute_progress_measures(Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction) override;
   [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double step_length) override;
   void postprocess_accepted_iterate(Iterate& iterate) override;

   [[nodiscard]] size_t get_hessian_evaluation_count() const override;
   [[nodiscard]] size_t get_number_subproblems_solved() const override;

private:
   const OptimalityProblem optimality_problem;
   l1RelaxedProblem feasibility_problem;
   std::unique_ptr<Subproblem> subproblem;
   const std::unique_ptr<GlobalizationStrategy> restoration_phase_strategy;
   const std::unique_ptr<GlobalizationStrategy> optimality_phase_strategy;
   Phase current_phase{Phase::OPTIMALITY};
   // statistics table
   int statistics_restoration_phase_column_order;

   [[nodiscard]] const NonlinearProblem& current_reformulated_problem() const;
   [[nodiscard]] GlobalizationStrategy& get_current_globalization_strategy() const;
   [[nodiscard]] Direction solve_optimality_problem(Statistics& statistics, Iterate& current_iterate);
   void switch_to_feasibility_restoration(Iterate& current_iterate);
   void switch_to_optimality(Iterate& current_iterate, Iterate& trial_iterate);

   // progress measures and their local models
   void set_infeasibility_measure(Iterate& iterate);
   [[nodiscard]] double generate_predicted_infeasibility_reduction_model(const Iterate& current_iterate,
         const Direction& direction, double step_length) const;
   void set_scaled_optimality_measure(Iterate& iterate);
   [[nodiscard]] std::function<double (double)> generate_predicted_scaled_optimality_reduction_model(const Iterate& current_iterate,
         const Direction& direction, double step_length) const;
};

#endif //UNO_FEASIBILITYRESTORATION_H