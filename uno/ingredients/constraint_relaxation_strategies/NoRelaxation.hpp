// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_NORELAXATION_H
#define UNO_NORELAXATION_H

#include <memory>
#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/globalization_strategies/MeritFunction.hpp"
#include "optimization/OptimizationProblem.hpp"

namespace uno {
   // forward declarations
   class Direction;
   class InequalityHandlingMethod;

   class NoRelaxation : public ConstraintRelaxationStrategy {
   public:
      NoRelaxation(const Model& model, Options& options);
      ~NoRelaxation() override;

      void initialize(Statistics& statistics, Iterate& initial_iterate, bool uses_trust_region,
         EvaluationCache& evaluation_cache, Options& options) override;

      // direction computation
      const Direction& compute_direction(Statistics& statistics, Iterate& current_iterate, double trust_region_radius,
         Evaluations& current_evaluations, WarmstartInformation& warmstart_information) override;
      [[nodiscard]] bool solving_feasibility_problem() const override;
      [[nodiscard]] bool test_infeasible_stationarity(Iterate& current_iterate, Evaluations& current_evaluations) const override;
      void switch_to_feasibility_problem(Statistics& statistics, Iterate& current_iterate, Evaluations& current_evaluations,
         WarmstartInformation& warmstart_information) override;

      // second-order corrections
      [[nodiscard]] bool has_second_order_corrections() const override;
      void initialize_second_order_corrections(const Iterate& current_iterate, const Iterate& trial_iterate,
         Evaluations& current_evaluations, Evaluations& trial_evaluations) override;
      const Direction& compute_second_order_correction(Iterate& current_iterate) override;
      void update_second_order_corrections(const Iterate& trial_iterate, Evaluations& trial_evaluations) override;

      // trial iterate acceptance
      [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Iterate& trial_iterate, const Direction& direction, double step_length, bool uses_trust_region,
         Evaluations& current_evaluations, Evaluations& trial_evaluations, WarmstartInformation& warmstart_information,
         UserCallbacks& user_callbacks) override;

      [[nodiscard]] std::string get_name() const override;

   private:
      const OptimizationProblem original_problem;
      std::unique_ptr<InequalityHandlingMethod> inequality_handling_method;
      MeritFunction globalization_strategy;
      Vector<double> initial_point;
   };
} // namespace

#endif //UNO_NORELAXATION_H