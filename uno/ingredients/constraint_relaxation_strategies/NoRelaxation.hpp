// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_NORELAXATION_H
#define UNO_NORELAXATION_H

#include <memory>
#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/globalization_strategies/l1MeritFunction.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethod.hpp"
#include "optimization/OptimizationProblem.hpp"

namespace uno {
   class NoRelaxation : public ConstraintRelaxationStrategy {
   public:
      NoRelaxation(const Model& model, const Options& options);
      ~NoRelaxation() override = default;

      void initialize(Statistics& statistics, const Model& model, Iterate& initial_iterate, Direction& direction,
         double trust_region_radius) override;

      // direction computation
      void compute_feasible_direction(Statistics& statistics, Iterate& current_iterate, Direction& direction,
         double trust_region_radius, WarmstartInformation& warmstart_information) override;
      [[nodiscard]] bool solving_feasibility_problem() const override;
      void switch_to_feasibility_problem(Statistics& statistics, Iterate& current_iterate, double trust_region_radius,
         WarmstartInformation& warmstart_information) override;

      // trial iterate acceptance
      [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Iterate& trial_iterate, const Direction& direction, double step_length, WarmstartInformation& warmstart_information,
         UserCallbacks& user_callbacks) override;
      [[nodiscard]] SolutionStatus check_termination(const Model& model, Iterate& iterate) override;

      [[nodiscard]] std::string get_name() const override;
      [[nodiscard]] size_t get_number_subproblems_solved() const override;

   private:
      const OptimizationProblem problem;
      std::unique_ptr<InequalityHandlingMethod> inequality_handling_method;
      std::unique_ptr<HessianModel> hessian_model;
      std::unique_ptr<InertiaCorrectionStrategy<double>> inertia_correction_strategy;
      l1MeritFunction globalization_strategy;
   };
} // namespace

#endif //UNO_NORELAXATION_H