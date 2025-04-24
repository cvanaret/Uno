// Copyright (c) 2018-2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FEASIBILITYRESTORATION_H
#define UNO_FEASIBILITYRESTORATION_H

#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/globalization_strategies/ProgressMeasures.hpp"
#include "layers/SubproblemLayer.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   // forward declaration
   class OptimizationProblem;

   enum class Phase {FEASIBILITY_RESTORATION = 1, OPTIMALITY = 2};

   class FeasibilityRestoration : public ConstraintRelaxationStrategy {
   public:
      FeasibilityRestoration(size_t number_bound_constraints, const Options& options);
      ~FeasibilityRestoration() override = default;

      void initialize(Statistics& statistics, const Model& model, Iterate& initial_iterate, Direction& direction,
         const Options& options) override;

      // direction computation
      void compute_feasible_direction(Statistics& statistics, GlobalizationStrategy& globalization_strategy, const Model& model,
         Iterate& current_iterate, Direction& direction, double trust_region_radius, WarmstartInformation& warmstart_information) override;
      [[nodiscard]] bool solving_feasibility_problem() const override;
      void switch_to_feasibility_problem(Statistics& statistics, GlobalizationStrategy& globalization_strategy, const Model& model,
         Iterate& current_iterate, WarmstartInformation& warmstart_information) override;

      // trial iterate acceptance
      [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const Model& model, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, double step_length,
         WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) override;

      // primal-dual residuals
      void compute_primal_dual_residuals(const Model& model, Iterate& iterate) override;
      void set_dual_residuals_statistics(Statistics& statistics, const Iterate& iterate) const override;

      [[nodiscard]] std::string get_name() const override;
      [[nodiscard]] size_t get_hessian_evaluation_count() const override;
      [[nodiscard]] size_t get_number_subproblems_solved() const override;

   private:
      Phase current_phase{Phase::OPTIMALITY};
      const double constraint_violation_coefficient;
      SubproblemLayer optimality_subproblem_layer;
      SubproblemLayer feasibility_subproblem_layer;
      std::unique_ptr<InequalityHandlingMethod> optimality_inequality_handling_method;
      std::unique_ptr<InequalityHandlingMethod> feasibility_inequality_handling_method;
      const double linear_feasibility_tolerance;
      const bool switch_to_optimality_requires_linearized_feasibility;
      ProgressMeasures reference_optimality_progress{};
      Vector<double> reference_optimality_primals{};

      void solve_subproblem(Statistics& statistics, InequalityHandlingMethod& inequality_handling_method, const OptimizationProblem& problem,
         Iterate& current_iterate, const Multipliers& current_multipliers, Direction& direction, SubproblemLayer& subproblem_layer,
         double trust_region_radius, WarmstartInformation& warmstart_information);
      void switch_to_optimality_phase(Iterate& current_iterate, GlobalizationStrategy& globalization_strategy, const Model& model,
         Iterate& trial_iterate, WarmstartInformation& warmstart_information);

      void evaluate_progress_measures(const OptimizationProblem& problem, InequalityHandlingMethod& inequality_handling_method,
         const Model& model, Iterate& iterate) const override;
      [[nodiscard]] ProgressMeasures compute_predicted_reductions(InequalityHandlingMethod& inequality_handling_method,
         const Model& model, const Iterate& current_iterate, const Direction& direction, double step_length) const;
      [[nodiscard]] bool can_switch_to_optimality_phase(const Iterate& current_iterate, const GlobalizationStrategy& globalization_strategy,
         const Model& model, const Iterate& trial_iterate, const Direction& direction, double step_length);
   };
} // namespace

#endif //UNO_FEASIBILITYRESTORATION_H