// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_UNCONSTRAINEDSTRATEGY_H
#define UNO_UNCONSTRAINEDSTRATEGY_H

#include <memory>
#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/globalization_strategies/ProgressMeasures.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethod.hpp"
#include "layers/SubproblemLayer.hpp"

namespace uno {
   // forward declarations
   class OptimizationProblem;

   class UnconstrainedStrategy : public ConstraintRelaxationStrategy {
   public:
      UnconstrainedStrategy(size_t number_bound_constraints, const Options& options);
      ~UnconstrainedStrategy() override = default;

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

   private:
      std::unique_ptr<InequalityHandlingMethod> inequality_handling_method;
      SubproblemLayer subproblem_layer;

      void solve_subproblem(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,
         const Multipliers& current_multipliers, Direction& direction, double trust_region_radius,
         WarmstartInformation& warmstart_information);

      void evaluate_progress_measures(const OptimizationProblem& problem, InequalityHandlingMethod& inequality_handling_method,
         const Model& model, Iterate& iterate) const override;
      [[nodiscard]] ProgressMeasures compute_predicted_reduction_models(InequalityHandlingMethod& inequality_handling_method,
         const Model& model, const Iterate& current_iterate, const Direction& direction, double step_length) const;
   };
} // namespace

#endif //UNO_UNCONSTRAINEDSTRATEGY_H
