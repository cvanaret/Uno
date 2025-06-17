// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_L1RELAXATION_H
#define UNO_L1RELAXATION_H

#include <memory>
#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/regularization_strategies/RegularizationStrategy.hpp"
#include "optimization/Multipliers.hpp"

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
      l1Relaxation(size_t number_bound_constraints, const Options& options);
      ~l1Relaxation() override = default;

      void initialize(Statistics& statistics, const Model& model, Iterate& initial_iterate, Direction& direction,
         const Options& options) override;

      // direction computation
      void compute_feasible_direction(Statistics& statistics, GlobalizationStrategy& globalization_strategy, const Model& model,
         Iterate& current_iterate, Direction& direction, double trust_region_radius, WarmstartInformation& warmstart_information) override;
      [[nodiscard]] bool solving_feasibility_problem() const override;
      void switch_to_feasibility_problem(Statistics& statistics, GlobalizationStrategy& globalization_strategy, const Model& model,
         Iterate& current_iterate, WarmstartInformation& warmstart_information) override;

      // trial iterate acceptance
      [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy, const Model& model,
         Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, double step_length,
         WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) override;

      // primal-dual residuals
      void compute_primal_dual_residuals(const Model& model, Iterate& iterate) override;
      void set_dual_residuals_statistics(Statistics& statistics, const Iterate& iterate) const override;

      [[nodiscard]] std::string get_name() const override;
      [[nodiscard]] size_t get_hessian_evaluation_count() const override;
      [[nodiscard]] size_t get_number_subproblems_solved() const override;

   protected:
      double penalty_parameter;
      const double constraint_violation_coefficient;
      std::unique_ptr<HessianModel> l1_relaxed_hessian_model;
      std::unique_ptr<HessianModel> feasibility_hessian_model;
      std::unique_ptr<RegularizationStrategy<double>> l1_relaxed_regularization_strategy;
      std::unique_ptr<RegularizationStrategy<double>> feasibility_regularization_strategy;
      std::unique_ptr<InequalityHandlingMethod> inequality_handling_method;
      std::unique_ptr<InequalityHandlingMethod> feasibility_inequality_handling_method;
      const double tolerance;
      const l1RelaxationParameters parameters;
      const double small_duals_threshold;
      // preallocated temporary multipliers
      Multipliers trial_multipliers{};

      void solve_sequence_of_relaxed_subproblems(Statistics& statistics, const Model& model, Iterate& current_iterate, Direction& direction,
         double trust_region_radius, WarmstartInformation& warmstart_information);
      void solve_l1_relaxed_problem(Statistics& statistics, const Model& model, Iterate& current_iterate, Direction& direction,
         double current_penalty_parameter, double trust_region_radius, WarmstartInformation& warmstart_information);
      void solve_subproblem(Statistics& statistics, InequalityHandlingMethod& inequality_handling_method, const OptimizationProblem& problem,
         Iterate& current_iterate, const Multipliers& current_multipliers, Direction& direction, HessianModel& hessian_model,
         RegularizationStrategy<double>& regularization_strategy, double trust_region_radius, WarmstartInformation& warmstart_information);

      // functions that decrease the penalty parameter to enforce particular conditions
      void decrease_parameter_aggressively(const Model& model, Iterate& current_iterate, const Direction& direction);
      double compute_infeasible_dual_error(const Model& model, Iterate& current_iterate) const;
      void enforce_linearized_residual_sufficient_decrease(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Direction& direction, double linearized_residual, double residual_lowest_violation, double trust_region_radius,
         WarmstartInformation& warmstart_information);
      [[nodiscard]] bool linearized_residual_sufficient_decrease(const Iterate& current_iterate, double linearized_residual,
         double residual_lowest_violation) const;
      void enforce_descent_direction_for_l1_merit(Statistics& statistics, const Model& model, Iterate& current_iterate, Direction& direction,
         const Direction& feasibility_direction, double trust_region_radius, WarmstartInformation& warmstart_information);
      [[nodiscard]] bool is_descent_direction_for_l1_merit_function(const Iterate& current_iterate, const Direction& direction,
         const Direction& feasibility_direction) const;

      void evaluate_progress_measures(InequalityHandlingMethod& inequality_handling_method, const Model& model, Iterate& iterate) const override;

      void check_exact_relaxation(Iterate& iterate) const;
   };
} // namespace

#endif //UNO_L1RELAXATION_H
