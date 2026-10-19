// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INEQUALITYHANDLINGMETHOD_H
#define UNO_INEQUALITYHANDLINGMETHOD_H

#include <optional>
#include <string>
#include "ingredients/globalization_strategies/PredictedReductionModels.hpp"
#include "linear_algebra/Norm.hpp"

namespace uno {
   // forward declarations
   class Direction;
   class Evaluations;
   class GlobalizationStrategy;
   class Iterate;
   class l1RelaxedProblem;
   class Multipliers;
   class OptimizationProblem;
   class Options;
   class ProgressMeasures;
   class Statistics;
   class SolverWorkspace;
   class Subproblem;
   template <typename T>
   class Vector;
   class WarmstartInformation;

   // the PhaseHandoff class is an object that is exchanged between phases. It is exported from the optimality
   // InequalityHandlingMethod and imported by the restoration InequalityHandlingMethod in FeasibilityRestoration
   struct PhaseHandoff {
      std::optional<double> barrier_parameter; // set by barrier methods only
   };
   
   class InequalityHandlingMethod {
   public:
      InequalityHandlingMethod(const OptimizationProblem& problem, const Options& options);
      virtual ~InequalityHandlingMethod() = default;

      virtual void initialize_memory() = 0;
      virtual void initialize_statistics(Statistics& statistics) = 0;
      virtual void create_initial_iterate(Iterate& iterate, Evaluations& evaluations, double multipliers_threshold) const = 0;
      [[nodiscard]] virtual std::pair<size_t, size_t> get_dimensions() const = 0;
      virtual void postprocess_solution(Iterate& iterate, Evaluations& evaluations) const = 0;

      [[nodiscard]] virtual bool update_parameterization(Statistics& statistics, const Iterate& current_iterate,
         Evaluations& evaluations) = 0;
      // unconditional update of the parameterization (e.g. after tiny steps). Returns true if it changed
      [[nodiscard]] virtual bool force_parameterization_update(Statistics& statistics) = 0;
      [[nodiscard]] virtual const Direction& solve(Statistics& statistics, const Iterate& current_iterate,
         double trust_region_radius, const Vector<double>& initial_point, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) = 0;

      virtual void import_handoff(const PhaseHandoff& handoff, Iterate& current_iterate) = 0;
      [[nodiscard]] virtual PhaseHandoff export_handoff() const = 0;
      // upon a phase switch, the incoming iterate may lie outside the (dynamically relaxed) bounds of this phase:
      // methods that maintain internal bounds relax them to strictly contain the iterate. No-op by default
      virtual void relax_bounds_around(const Iterate& /*iterate*/) const { }
      virtual void set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate, Evaluations& evaluations) = 0;
      [[nodiscard]] virtual double proximal_coefficient() const = 0;

      [[nodiscard]] virtual bool has_second_order_corrections() const = 0;
      virtual void initialize_second_order_corrections(const Iterate& current_iterate, const Iterate& trial_iterate,
         Evaluations& current_evaluations, Evaluations& trial_evaluations) = 0;
      [[nodiscard]] virtual const Direction& compute_second_order_correction(const Iterate& current_iterate) = 0;
      virtual void update_second_order_corrections(const Iterate& trial_iterate, Evaluations& trial_evaluations) = 0;

      virtual void compute_bound_dual_direction(const Vector<double>& current_primals, const Multipliers& current_multipliers,
         const Vector<double>& direction_primals, Multipliers& direction_multipliers, double& bound_dual_step_length) const = 0;

      virtual void evaluate_progress_measures(Iterate& iterate, Evaluations& evaluations) const = 0;
      [[nodiscard]] virtual PredictedReductionModels build_predicted_reduction_models(const Iterate& current_iterate,
         const Direction& direction, Evaluations& current_evaluations) const = 0;
      [[nodiscard]] virtual bool is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, Evaluations& trial_evaluations,
         const ProgressMeasures& predicted_reductions) const = 0;
      // accept the trial iterate without globalization test (tiny steps)
      virtual void accept_iterate_unconditionally(Iterate& trial_iterate, Evaluations& trial_evaluations) const = 0;
      virtual void notify_trial_iterate(Statistics& statistics, const Iterate& current_iterate, const Iterate& trial_iterate,
         Evaluations& current_evaluations, Evaluations& trial_evaluations) = 0;

      virtual void compute_residuals(Iterate& iterate, Evaluations& evaluations) const = 0;
      // square problem (as many free variables as equality constraints, IPOPT's IsSquareProblem): the constraint
      // multipliers are recomputed as least-squares estimates at each accepted iterate (IPOPT's ComputeFeasibilityMultipliers)
      virtual void update_multipliers_of_square_problem(Iterate& /*iterate*/, Evaluations& /*evaluations*/,
         double /*multipliers_threshold*/) { }
      // primal infeasibility of the (reformulated) problem, as used in the termination tests
      [[nodiscard]] virtual double compute_primal_infeasibility(const Iterate& iterate, Evaluations& evaluations) const = 0;

      [[nodiscard]] virtual std::string get_name() const = 0;

   protected:
      const OptimizationProblem& problem;
      const Norm progress_norm;
      const Norm residual_norm;

      void evaluate_progress_measures(const OptimizationProblem& problem, Iterate& iterate, Evaluations& evaluations) const;
      void compute_residuals(const OptimizationProblem& problem, Iterate& iterate, Evaluations& evaluations) const;
      bool is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const Subproblem& subproblem, const Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         Evaluations& trial_evaluations, const ProgressMeasures& predicted_reductions) const;
      void accept_iterate_unconditionally(const Subproblem& subproblem, Iterate& trial_iterate, Evaluations& trial_evaluations) const;
   };
} // namespace

#endif // UNO_INEQUALITYHANDLINGMETHOD_H