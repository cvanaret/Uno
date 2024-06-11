// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FEASIBILITYRESTORATION_H
#define UNO_FEASIBILITYRESTORATION_H

#include <memory>
#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/globalization_strategy/GlobalizationStrategy.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "reformulation/OptimalityProblem.hpp"
#include "reformulation/l1RelaxedProblem.hpp"

enum class Phase {FEASIBILITY_RESTORATION = 1, OPTIMALITY = 2};

class FeasibilityRestoration : public ConstraintRelaxationStrategy {
public:
   FeasibilityRestoration(const Model& model, const Options& options);

   void initialize(Statistics& statistics, Iterate& initial_iterate, const Options& options) override;
   void set_trust_region_radius(double trust_region_radius) override;

   [[nodiscard]] size_t maximum_number_variables() const override;
   [[nodiscard]] size_t maximum_number_constraints() const override;

   // direction computation
   void compute_feasible_direction(Statistics& statistics, Iterate& current_iterate, Direction& direction, WarmstartInformation& warmstart_information) override;
   void compute_feasible_direction(Statistics& statistics, Iterate& current_iterate, Direction& direction, const Vector<double>& initial_point,
         WarmstartInformation& warmstart_information) override;
   [[nodiscard]] bool solving_feasibility_problem() const override;
   void switch_to_feasibility_problem(Statistics& statistics, Iterate& current_iterate) override;

   // trial iterate acceptance
   void compute_progress_measures(Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, double step_length) override;
   [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double step_length) override;

   // primal-dual residuals
   void compute_primal_dual_residuals(Iterate& iterate) override;

   [[nodiscard]] size_t get_hessian_evaluation_count() const override;
   [[nodiscard]] size_t get_number_subproblems_solved() const override;

private:
   const OptimalityProblem optimality_problem;
   l1RelaxedProblem feasibility_problem;
   std::unique_ptr<Subproblem> subproblem;
   const std::unique_ptr<GlobalizationStrategy> globalization_strategy;
   Phase current_phase{Phase::OPTIMALITY};
   const double linear_feasibility_tolerance;
   const bool switch_to_optimality_requires_acceptance;
   const bool switch_to_optimality_requires_linearized_feasibility;
   bool switching_to_optimality_phase{false};

   [[nodiscard]] const OptimizationProblem& current_problem() const;
   void solve_subproblem(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate, const Multipliers& current_multipliers,
         Direction& direction, WarmstartInformation& warmstart_information);
   void switch_to_optimality_phase(Iterate& current_iterate, Iterate& trial_iterate);

   void evaluate_progress_measures(Iterate& iterate) const;
   [[nodiscard]] ProgressMeasures compute_predicted_reduction_models(Iterate& current_iterate, const Direction& direction, double step_length);

   void set_dual_residuals_statistics(Statistics& statistics, const Iterate& iterate) const override;
};

#endif //UNO_FEASIBILITYRESTORATION_H
