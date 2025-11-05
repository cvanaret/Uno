// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_L1RELAXATION_H
#define UNO_L1RELAXATION_H

#include <memory>
#include "ConstraintRelaxationStrategy.hpp"
#include "l1RelaxedProblem.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethod.hpp"

namespace uno {
   struct l1RelaxationParameters {
      double beta;
      double theta;
      double kappa_rho;
      double kappa_lambda;
      double epsilon;
      double omega;
      double delta;
   };

   class l1Relaxation : public ConstraintRelaxationStrategy {
   public:
      explicit l1Relaxation(const Options& options);
      ~l1Relaxation() override = default;

      void initialize(Statistics& statistics, const Model& model, Iterate& initial_iterate, Direction& direction,
         double trust_region_radius, const Options& options) override;

      // direction computation
      void compute_feasible_direction(Statistics& statistics, GlobalizationStrategy& globalization_strategy, Iterate& current_iterate,
         Direction& direction, double trust_region_radius, WarmstartInformation& warmstart_information) override;
      [[nodiscard]] bool solving_feasibility_problem() const override;
      void switch_to_feasibility_problem(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         Iterate& current_iterate, double trust_region_radius, WarmstartInformation& warmstart_information) override;

      // trial iterate acceptance
      [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         double trust_region_radius, const Model& model, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         double step_length, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) override;
      [[nodiscard]] SolutionStatus check_termination(const Model& model, Iterate& iterate) override;

      [[nodiscard]] std::string get_name() const override;
      [[nodiscard]] size_t get_number_subproblems_solved() const override;

   protected:
      std::unique_ptr<l1RelaxedProblem> relaxed_problem{};
      std::unique_ptr<const l1RelaxedProblem> feasibility_problem{};
      std::unique_ptr<InequalityHandlingMethod> inequality_handling_method;
      std::unique_ptr<InequalityHandlingMethod> feasibility_inequality_handling_method;
      std::unique_ptr<HessianModel> hessian_model{};
      std::unique_ptr<InertiaCorrectionStrategy<double>> inertia_correction_strategy;
      Multipliers feasibility_multipliers;
      double penalty_parameter;
      const double tolerance;
      const l1RelaxationParameters parameters;

      // auxiliary functions with the notations of the paper
      [[nodiscard]] double v(const Iterate& current_iterate) const;
      [[nodiscard]] double l(const Direction& direction, const Iterate& current_iterate) const;
      [[nodiscard]] double delta_l(const Direction& direction, const Iterate& current_iterate) const;
      [[nodiscard]] double Ropt(Iterate& current_iterate, double objective_multiplier, const Multipliers& multipliers) const;
      [[nodiscard]] double Rinf(Iterate& current_iterate, const Multipliers& multipliers) const;
      [[nodiscard]] double delta_m(const Direction& direction, const Iterate& current_iterate, double objective_multiplier) const;
      [[nodiscard]] double compute_zeta(const Direction& direction, const Iterate& current_iterate) const;
      [[nodiscard]] double compute_w(const Direction& feasibility_direction, const Direction& optimality_direction, const Iterate& current_iterate);
      void check_exact_relaxation(const Iterate& iterate) const;
   };
} // namespace

#endif //UNO_L1RELAXATION_H