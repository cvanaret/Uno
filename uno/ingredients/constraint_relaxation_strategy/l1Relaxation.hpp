// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_L1RELAXATION_H
#define UNO_L1RELAXATION_H

#include <memory>
#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/globalization_strategy/GlobalizationStrategy.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "reformulation/l1RelaxedProblem.hpp"

struct l1RelaxationParameters {
   bool fixed_parameter;
   double decrease_factor;
   double epsilon1;
   double epsilon2;
   double residual_small_threshold;
};

class l1Relaxation : public ConstraintRelaxationStrategy {
public:
   l1Relaxation(const Model& model, const Options& options);

   void initialize(Statistics& statistics, Iterate& initial_iterate, const Options& options) override;
   void set_trust_region_radius(double trust_region_radius) override;

   [[nodiscard]] size_t maximum_number_variables() const override;
   [[nodiscard]] size_t maximum_number_constraints() const override;

   // direction computation
   void compute_feasible_direction(Statistics& statistics, Iterate& current_iterate, Direction& direction,
         WarmstartInformation& warmstart_information) override;
   void compute_feasible_direction(Statistics& statistics, Iterate& current_iterate, Direction& direction, const std::vector<double>& initial_point,
         WarmstartInformation& warmstart_information) override;
   [[nodiscard]] bool solving_feasibility_problem() const override;
   void switch_to_feasibility_problem(Statistics& statistics, Iterate& current_iterate) override;

   // trial iterate acceptance
   void compute_progress_measures(Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, double step_length) override;
   [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double step_length) override;

   [[nodiscard]] size_t get_hessian_evaluation_count() const override;
   [[nodiscard]] size_t get_number_subproblems_solved() const override;

protected:
   const l1RelaxedProblem feasibility_problem;
   l1RelaxedProblem l1_relaxed_problem;
   std::unique_ptr<Subproblem> subproblem;
   const std::unique_ptr<GlobalizationStrategy> globalization_strategy;
   double penalty_parameter;
   const double tolerance;
   const l1RelaxationParameters parameters;
   const double small_duals_threshold;
   // preallocated temporary multipliers
   Multipliers trial_multipliers;

   void solve_sequence_of_relaxed_subproblems(Statistics& statistics, Iterate& current_iterate, Direction& direction,
         WarmstartInformation& warmstart_information);
   void solve_l1_relaxed_problem(Statistics& statistics, Iterate& current_iterate, Direction& direction, double current_penalty_parameter,
         const WarmstartInformation& warmstart_information);
   void solve_subproblem(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate, Direction& direction,
         const WarmstartInformation& warmstart_information);

   // functions that decrease the penalty parameter to enforce particular conditions
   void decrease_parameter_aggressively(Iterate& current_iterate, const Direction& direction);
   double compute_infeasible_dual_error(Iterate& current_iterate);
   void enforce_linearized_residual_sufficient_decrease(Statistics& statistics, Iterate& current_iterate, Direction& direction,
         double linearized_residual, double residual_lowest_violation, WarmstartInformation& warmstart_information);
   [[nodiscard]] bool linearized_residual_sufficient_decrease(const Iterate& current_iterate, double linearized_residual,
         double residual_lowest_violation) const;
   void enforce_descent_direction_for_l1_merit(Statistics& statistics, Iterate& current_iterate, Direction& direction,
         const Direction& direction_lowest_violation, WarmstartInformation& warmstart_information);
   [[nodiscard]] bool is_descent_direction_for_l1_merit_function(const Iterate& current_iterate, const Direction& direction,
         const Direction& direction_lowest_violation) const;

   void evaluate_progress_measures(Iterate& iterate) const;
   [[nodiscard]] ProgressMeasures compute_predicted_reduction_models(Iterate& current_iterate, const Direction& direction, double step_length);

   void set_dual_residuals_statistics(Statistics& statistics, const Iterate& iterate) const override;
   void check_exact_relaxation(Iterate& iterate) const;
};

#endif //UNO_L1RELAXATION_H
