// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_L1RELAXATION_H
#define UNO_L1RELAXATION_H

#include <memory>
#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/globalization_strategies/ProgressMeasures.hpp"
#include "optimization/Multipliers.hpp"
#include "l1RelaxedProblem.hpp"

namespace uno {
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

      [[nodiscard]] size_t maximum_number_variables() const override;
      [[nodiscard]] size_t maximum_number_constraints() const override;

      // direction computation
      void compute_feasible_direction(Statistics& statistics, Iterate& current_iterate, Direction& direction, double trust_region_radius,
            WarmstartInformation& warmstart_information) override;
      [[nodiscard]] bool solving_feasibility_problem() const override;
      void switch_to_feasibility_problem(Statistics& statistics, Iterate& current_iterate, WarmstartInformation& warmstart_information) override;

      // trial iterate acceptance
      [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
            double step_length, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) override;

      // primal-dual residuals
      void compute_primal_dual_residuals(Iterate& iterate) override;
      void set_dual_residuals_statistics(Statistics& statistics, const Iterate& iterate) const override;

   protected:
      const l1RelaxedProblem feasibility_problem;
      l1RelaxedProblem l1_relaxed_problem;
      double penalty_parameter;
      const double tolerance;
      const l1RelaxationParameters parameters;
      const double small_duals_threshold;
      // preallocated temporary multipliers
      Multipliers trial_multipliers;

      // delegating constructor
      l1Relaxation(const Model& model, l1RelaxedProblem&& feasibility_problem, l1RelaxedProblem&& l1_relaxed_problem, const Options& options);

      void solve_sequence_of_relaxed_subproblems(Statistics& statistics, Iterate& current_iterate, Direction& direction, double trust_region_radius,
            WarmstartInformation& warmstart_information);
      void solve_l1_relaxed_problem(Statistics& statistics, Iterate& current_iterate, Direction& direction, double trust_region_radius,
         double current_penalty_parameter, WarmstartInformation& warmstart_information);
      void solve_subproblem(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate, const Multipliers& current_multipliers,
            Direction& direction, double trust_region_radius, WarmstartInformation& warmstart_information);

      // functions that decrease the penalty parameter to enforce particular conditions
      void decrease_parameter_aggressively(Iterate& current_iterate, const Direction& direction);
      double compute_infeasible_dual_error(Iterate& current_iterate);
      void enforce_linearized_residual_sufficient_decrease(Statistics& statistics, Iterate& current_iterate, Direction& direction,
            double trust_region_radius, double linearized_residual, double residual_lowest_violation, WarmstartInformation& warmstart_information);
      [[nodiscard]] bool linearized_residual_sufficient_decrease(const Iterate& current_iterate, double linearized_residual,
            double residual_lowest_violation) const;
      void enforce_descent_direction_for_l1_merit(Statistics& statistics, Iterate& current_iterate, Direction& direction,
            const Direction& feasibility_direction, double trust_region_radius, WarmstartInformation& warmstart_information);
      [[nodiscard]] bool is_descent_direction_for_l1_merit_function(const Iterate& current_iterate, const Direction& direction,
            const Direction& feasibility_direction) const;

      void evaluate_progress_measures(Iterate& iterate) const override;
      [[nodiscard]] ProgressMeasures compute_predicted_reduction_models(Iterate& current_iterate, const Direction& direction, double step_length);

      void check_exact_relaxation(Iterate& iterate) const;
   };
} // namespace

#endif //UNO_L1RELAXATION_H
