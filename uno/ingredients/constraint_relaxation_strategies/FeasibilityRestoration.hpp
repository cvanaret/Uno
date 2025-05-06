// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FEASIBILITYRESTORATION_H
#define UNO_FEASIBILITYRESTORATION_H

#include <memory>
#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/globalization_strategies/ProgressMeasures.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   // forward declarations
   class OptimizationProblem;

   enum class Phase {FEASIBILITY_RESTORATION = 1, OPTIMALITY = 2};

   class FeasibilityRestoration : public ConstraintRelaxationStrategy {
   public:
      explicit FeasibilityRestoration(const Options& options);
      ~FeasibilityRestoration() override = default;

      void initialize(Statistics& statistics, const Model& model, Iterate& initial_iterate, Direction& direction, const Options& options) override;

      // direction computation
      void compute_feasible_direction(Statistics& statistics, const Model& model, Iterate& current_iterate, Direction& direction,
         WarmstartInformation& warmstart_information) override;
      [[nodiscard]] bool solving_feasibility_problem() const override;
      void switch_to_feasibility_problem(Statistics& statistics, const Model& model, Iterate& current_iterate,
         WarmstartInformation& warmstart_information) override;

      // trial iterate acceptance
      [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, const Model& model, Iterate& current_iterate, Iterate& trial_iterate,
         const Direction& direction, double step_length, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) override;

      // primal-dual residuals
      void compute_primal_dual_residuals(const Model& model, Iterate& iterate) override;
      void set_dual_residuals_statistics(Statistics& statistics, const Iterate& iterate) const override;

      [[nodiscard]] size_t get_hessian_evaluation_count() const override;

   private:
      Phase current_phase{Phase::OPTIMALITY};
      const double constraint_violation_coefficient;
      const bool convexify;
      const std::unique_ptr<HessianModel> optimality_hessian_model;
      const std::unique_ptr<HessianModel> feasibility_hessian_model;
      const double linear_feasibility_tolerance;
      const bool switch_to_optimality_requires_linearized_feasibility;
      ProgressMeasures reference_optimality_progress{};
      Vector<double> reference_optimality_primals{};

      void solve_subproblem(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate, const Multipliers& current_multipliers,
         Direction& direction, HessianModel& hessian_model, WarmstartInformation& warmstart_information);
      void switch_to_optimality_phase(Iterate& current_iterate, const Model& model, Iterate& trial_iterate,
         WarmstartInformation& warmstart_information);

      void evaluate_progress_measures(const Model& model, Iterate& iterate) const override;
      [[nodiscard]] ProgressMeasures compute_predicted_reduction_models(const Model& model, const Iterate& current_iterate,
         const Direction& direction, double step_length) const;
      [[nodiscard]] bool can_switch_to_optimality_phase(const Iterate& current_iterate, const Model& model, const Iterate& trial_iterate,
         const Direction& direction, double step_length);
   };
} // namespace

#endif //UNO_FEASIBILITYRESTORATION_H
