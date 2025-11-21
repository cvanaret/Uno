// Copyright (c) 2018-2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FEASIBILITYRESTORATION_H
#define UNO_FEASIBILITYRESTORATION_H

#include <memory>
#include "ConstraintRelaxationStrategy.hpp"
#include "relaxed_problems/l1RelaxedProblem.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/globalization_strategies/l1MeritFunction.hpp"
#include "ingredients/globalization_strategies/ProgressMeasures.hpp"
#include "ingredients/inertia_correction_strategies/InertiaCorrectionStrategy.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   // forward declaration
   class InequalityHandlingMethod;

   enum class Phase {FEASIBILITY_RESTORATION = 1, OPTIMALITY = 2};

   class FeasibilityRestoration : public ConstraintRelaxationStrategy {
   public:
      FeasibilityRestoration(const Model& model, bool use_trust_region, const Options& options);
      ~FeasibilityRestoration() override = default;

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
      Phase current_phase{Phase::OPTIMALITY};
      const double constraint_violation_coefficient;
      const OptimizationProblem optimality_problem;
      l1RelaxedProblem feasibility_problem;
      std::unique_ptr<HessianModel> optimality_hessian_model;
      std::unique_ptr<HessianModel> feasibility_hessian_model;
      std::unique_ptr<InertiaCorrectionStrategy<double>> optimality_inertia_correction_strategy;
      std::unique_ptr<InertiaCorrectionStrategy<double>> feasibility_inertia_correction_strategy;
      std::unique_ptr<InequalityHandlingMethod> optimality_inequality_handling_method;
      std::unique_ptr<InequalityHandlingMethod> feasibility_inequality_handling_method;
      std::unique_ptr<GlobalizationStrategy> optimality_globalization_strategy;
      l1MeritFunction feasibility_globalization_strategy;
      // the class maintains multipliers for the other phase (feasibility multipliers if we are in the optimality phase,
      // and vice versa). These multipliers and those of the iterate are swapped whenever we switch phases.
      Multipliers other_phase_multipliers;
      const double linear_feasibility_tolerance;
      const bool switch_to_optimality_requires_linearized_feasibility;
      ProgressMeasures reference_optimality_progress{};
      Vector<double> reference_optimality_primals{};
      bool first_switch_to_feasibility{true};

      void solve_subproblem(Statistics& statistics, InequalityHandlingMethod& inequality_handling_method, const OptimizationProblem& problem,
         Iterate& current_iterate, Direction& direction, double trust_region_radius, WarmstartInformation& warmstart_information);
      void switch_back_to_optimality_phase(Iterate& current_iterate, Iterate& trial_iterate);

      [[nodiscard]] bool can_switch_to_optimality_phase(const Iterate& current_iterate, const Model& model,
         const Iterate& trial_iterate, const Direction& direction, double step_length) const;
   };
} // namespace

#endif //UNO_FEASIBILITYRESTORATION_H