// Copyright (c) 2018-2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FEASIBILITYRESTORATION_H
#define UNO_FEASIBILITYRESTORATION_H

#include <memory>
#include "ConstraintRelaxationStrategy.hpp"
#include "relaxed_problems/l1RelaxedProblem.hpp"
#include "ingredients/globalization_strategies/MeritFunction.hpp"
#include "ingredients/globalization_strategies/ProgressMeasures.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   // forward declaration
   class Direction;
   class HessianModel;
   class InequalityHandlingMethod;
   class InertiaCorrectionStrategy;
   class SubproblemSolver;

   enum class Phase {FEASIBILITY_RESTORATION = 1, OPTIMALITY = 2};

   class FeasibilityRestoration : public ConstraintRelaxationStrategy {
   public:
      FeasibilityRestoration(const Model& model, bool use_trust_region, Options& options);
      ~FeasibilityRestoration() override;

      void initialize(Statistics& statistics, Iterate& initial_iterate, bool uses_trust_region,
         EvaluationCache& evaluation_cache, Options& options) override;

      // direction computation
      Direction& compute_feasible_direction(Statistics& statistics, Iterate& current_iterate, double trust_region_radius,
         Evaluations& current_evaluations, WarmstartInformation& warmstart_information) override;
      [[nodiscard]] bool solving_feasibility_problem() const override;
      void switch_to_feasibility_problem(Statistics& statistics, Iterate& current_iterate, Evaluations& current_evaluations,
         WarmstartInformation& warmstart_information) override;

      // second-order corrections
      [[nodiscard]] bool has_second_order_corrections() const override;
      const Direction& compute_second_order_correction(Iterate& current_iterate, const Vector<double>& constraints_SOC) override;

      // trial iterate acceptance
      [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Iterate& trial_iterate, const Direction& direction, double step_length, bool uses_trust_region,
         Evaluations& current_evaluations, Evaluations& trial_evaluations, WarmstartInformation& warmstart_information,
         UserCallbacks& user_callbacks) override;

      [[nodiscard]] std::string get_name() const override;

   private:
      Phase current_phase{Phase::OPTIMALITY};
      const double constraint_violation_coefficient;
      const OptimizationProblem original_problem;
      l1RelaxedProblem feasibility_problem;
      std::unique_ptr<InequalityHandlingMethod> inequality_handling_method;
      std::unique_ptr<InequalityHandlingMethod> feasibility_inequality_handling_method;
      std::unique_ptr<GlobalizationStrategy> globalization_strategy;
      MeritFunction feasibility_globalization_strategy;
      Vector<double> initial_point;

      // the class maintains multipliers for the other phase (feasibility multipliers if we are in the optimality phase,
      // and vice versa). These multipliers and those of the iterate are swapped whenever we switch phases.
      Multipliers other_phase_multipliers;
      const double linear_feasibility_tolerance;
      const bool switch_to_optimality_requires_linearized_feasibility;
      ProgressMeasures reference_optimality_progress{};
      Vector<double> reference_optimality_primals{};

      Direction& solve_subproblem(Statistics& statistics, InequalityHandlingMethod& inequality_handling_method,
         GlobalizationStrategy& globalization_strategy, Iterate& current_iterate, double trust_region_radius,
         Evaluations& current_evaluations, const WarmstartInformation& warmstart_information);
      void switch_back_to_optimality_phase(Iterate& current_iterate, Iterate& trial_iterate);

      [[nodiscard]] bool can_switch_to_optimality_phase(const Model& model, const Iterate& trial_iterate,
         const Direction& direction, double step_length, Evaluations& current_evaluations) const;
   };
} // namespace

#endif //UNO_FEASIBILITYRESTORATION_H