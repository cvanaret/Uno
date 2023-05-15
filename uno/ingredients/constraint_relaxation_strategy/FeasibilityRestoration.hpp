// Copyright (c) 2018-2023 Charlie Vanaret
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
   FeasibilityRestoration(Statistics& statistics, const Model& model, const Options& options);
   void initialize(Iterate& initial_iterate) override;

   void set_trust_region_radius(double trust_region_radius) override;

   // direction computation
   [[nodiscard]] Direction compute_feasible_direction(Statistics& statistics, Iterate& current_iterate,
         WarmstartInformation& warmstart_information) override;
   [[nodiscard]] Direction compute_feasible_direction(Statistics& statistics, Iterate& current_iterate, const std::vector<double>& initial_point,
         WarmstartInformation& warmstart_information) override;
   void switch_to_feasibility_problem(Iterate& current_iterate, WarmstartInformation& warmstart_information) override;

   // trial iterate acceptance
   void compute_progress_measures(Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, double step_length) override;
   [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double step_length) override;

   [[nodiscard]] size_t get_hessian_evaluation_count() const override;
   [[nodiscard]] size_t get_number_subproblems_solved() const override;

private:
   const OptimalityProblem optimality_problem;
   l1RelaxedProblem feasibility_problem;
   std::unique_ptr<Subproblem> subproblem;
   const std::unique_ptr<GlobalizationStrategy> restoration_phase_strategy;
   const std::unique_ptr<GlobalizationStrategy> optimality_phase_strategy;
   Phase current_phase{Phase::OPTIMALITY};
   const double l1_constraint_violation_coefficient;
   const double tolerance;
   const bool test_linearized_feasibility;
   bool switched_to_optimality_phase{false};

   [[nodiscard]] const NonlinearProblem& current_problem() const;
   [[nodiscard]] GlobalizationStrategy& current_globalization_strategy() const;
   [[nodiscard]] Direction solve_subproblem(Statistics& statistics, const NonlinearProblem& problem, Iterate& current_iterate,
         WarmstartInformation& warmstart_information);
   void switch_to_feasibility_restoration(Iterate& current_iterate, WarmstartInformation& warmstart_information);
   void switch_to_optimality(Iterate& current_iterate, Iterate& trial_iterate);

   // progress measures and their local models
   void set_progress_measures_for_optimality_problem(Iterate& iterate);
   void set_progress_measures_for_feasibility_problem(Iterate& iterate);
   [[nodiscard]] ProgressMeasures compute_predicted_reduction_models_for_optimality_problem(const Iterate& current_iterate, const Direction& direction,
         double step_length);
   [[nodiscard]] ProgressMeasures compute_predicted_reduction_models_for_feasibility_problem(const Iterate& current_iterate, const Direction& direction,
         double step_length);
};

#endif //UNO_FEASIBILITYRESTORATION_H