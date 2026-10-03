// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_NOINEQUALITYREFORMULATION_H
#define UNO_NOINEQUALITYREFORMULATION_H

#include <memory>
#include <string>
#include <vector>
#include "InequalityHandlingMethod.hpp"
#include "ingredients/hessian_models/HessianSubproblemSolverJointFactory.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "options/Options.hpp"

namespace uno {
   // forward declaration
   class HessianModel;
   class InertiaCorrectionStrategy;
   class SubproblemSolver;

   class NoInequalityReformulation : public InequalityHandlingMethod {
   public:
      NoInequalityReformulation(std::string name, const OptimizationProblem& problem, bool uses_trust_region,
         double objective_multiplier, const Options& options, std::vector<OptionOverride>& option_overrides);
      ~NoInequalityReformulation() override = default;

      void initialize_memory() override;
      void initialize_statistics(Statistics& statistics) override;
      void create_initial_iterate(Iterate& initial_iterate, Evaluations& evaluations, double multipliers_threshold) const override;
      [[nodiscard]] std::pair<size_t, size_t> get_dimensions() const override;

      [[nodiscard]] bool update_parameterization(Statistics& statistics, const Iterate& current_iterate,
         Evaluations& evaluations) override;
      [[nodiscard]] bool force_parameterization_update(Statistics& statistics) override;
      [[nodiscard]] const Direction& solve(Statistics& statistics, const Iterate& current_iterate, double trust_region_radius,
         const Vector<double>& initial_point, Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) override;

      void import_handoff(const PhaseHandoff& handoff, Iterate& current_iterate) override;
      [[nodiscard]] PhaseHandoff export_handoff() const override;
      void set_elastic_variable_values(const l1RelaxedProblem& feasibility_problem, Iterate& current_iterate,
         Evaluations& evaluations) override;
      [[nodiscard]] double proximal_coefficient() const override;

      [[nodiscard]] bool has_second_order_corrections() const override;
      void initialize_second_order_corrections(const Iterate& current_iterate, const Iterate& trial_iterate,
         Evaluations& current_evaluations, Evaluations& trial_evaluations) override;
      [[nodiscard]] const Direction& compute_second_order_correction(const Iterate& current_iterate) override;
      void update_second_order_corrections(const Iterate& trial_iterate, Evaluations& trial_evaluations) override;

      void compute_least_squares_multipliers(Iterate& iterate, Evaluations& evaluations) override;
      void compute_bound_dual_direction(const Vector<double>& current_primals, const Multipliers& current_multipliers,
         const Vector<double>& direction_primals, Multipliers& direction_multipliers, double& bound_dual_step_length) const override;

      void evaluate_progress_measures(Iterate& iterate, Evaluations& evaluations) const override;
      [[nodiscard]] PredictedReductionModels build_predicted_reduction_models(const Iterate& current_iterate,
         const Direction& direction, Evaluations& current_evaluations) const override;
      [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, Evaluations& trial_evaluations,
         const ProgressMeasures& predicted_reductions) const override;
      void accept_iterate_unconditionally(Iterate& trial_iterate, Evaluations& trial_evaluations) const override;
      void notify_trial_iterate(Statistics& statistics, const Iterate& current_iterate, const Iterate& trial_iterate,
         Evaluations& current_evaluations, Evaluations& trial_evaluations) override;

      void compute_residuals(Iterate& iterate, Evaluations& evaluations) const override;

      [[nodiscard]] std::string get_name() const override;

   protected:
      const std::string name;
      std::unique_ptr<InertiaCorrectionStrategy> inertia_correction_strategy;
      std::unique_ptr<HessianModel> hessian_model;
      std::unique_ptr<SubproblemSolver> subproblem_solver;
      std::unique_ptr<Subproblem> subproblem;
   };
} // namespace

#endif // UNO_NOINEQUALITYREFORMULATION_H