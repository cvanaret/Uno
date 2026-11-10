// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_CONSTRAINTRELAXATIONSTRATEGY_H
#define UNO_CONSTRAINTRELAXATIONSTRATEGY_H

#include <cstddef>
#include <optional>
#include <stdexcept>
#include <vector>
#include "ingredients/globalization_strategies/Decision.hpp"
#include "ingredients/globalization_strategies/PredictedReductionModels.hpp"
#include "linear_algebra/Norm.hpp"
#include "model/Model.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/SolutionStatus.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"

namespace uno {
   // forward declarations
   class Direction;
   class EvaluationCache;
   class GlobalizationStrategy;
   class OptimizationProblem;
   class Statistics;
   class UserCallbacks;
   class WarmstartInformation;

   class ConstraintRelaxationStrategy {
   public:
      explicit ConstraintRelaxationStrategy(const Options& options);
      virtual ~ConstraintRelaxationStrategy();

      virtual void initialize(Statistics& statistics, Iterate& initial_iterate, bool uses_trust_region,
         EvaluationCache& evaluation_cache, const Options& options, std::vector<OptionOverride>& option_overrides) = 0;
      virtual void postprocess_solution(Iterate& iterate, Evaluations& evaluations) const = 0;

      // direction computation
      virtual const Direction& compute_direction(Statistics& statistics, Iterate& current_iterate, double trust_region_radius,
         Evaluations& current_evaluations, WarmstartInformation& warmstart_information) = 0;
      [[nodiscard]] virtual bool solving_feasibility_problem() const = 0;
      virtual void switch_to_feasibility_problem(Statistics& statistics, Iterate& current_iterate, Evaluations& current_evaluations,
         WarmstartInformation& warmstart_information) = 0;

      // second-order corrections
      [[nodiscard]] virtual bool has_second_order_corrections() const = 0;
      virtual void initialize_second_order_corrections(const Iterate& current_iterate, const Iterate& trial_iterate,
         Evaluations& current_evaluations, Evaluations& trial_evaluations) = 0;
      virtual const Direction& compute_second_order_correction(Statistics& statistics, Iterate& current_iterate) = 0;
      virtual void update_second_order_corrections(const Iterate& trial_iterate, Evaluations& trial_evaluations) = 0;

      // trial iterate acceptance
      [[nodiscard]] virtual PredictedReductionModels build_predicted_reduction_models(const Iterate& current_iterate,
         const Direction& direction, Evaluations& current_evaluations) const = 0;
      [[nodiscard]] virtual Decision is_iterate_acceptable(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Iterate& trial_iterate, const Direction& direction, double step_length, bool uses_trust_region, Evaluations& current_evaluations,
         Evaluations& trial_evaluations, const ProgressMeasures& predicted_reductions, WarmstartInformation& warmstart_information,
         UserCallbacks& user_callbacks) = 0;
      // forwarded to the globalization strategy that tested the trial iterates of the step
      virtual void notify_step_limited_by(RejectionReason sole_rejection_reason) = 0;
      // accept a trial step without globalization test (IPOPT, Section 3.9), possibly after forcing an update of the
      // parameterization (e.g. a decrease of the barrier parameter). Returns true if the parameterization was updated
      [[nodiscard]] virtual bool accept_iterate_unconditionally(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Iterate& trial_iterate, const Direction& direction, double step_length, bool force_parameterization_update,
         Evaluations& current_evaluations, Evaluations& trial_evaluations, WarmstartInformation& warmstart_information,
         UserCallbacks& user_callbacks) = 0;
      // the line search failed. Returns true if the trial iterate was set to a terminal iterate (the current iterate or
      // the last acceptable iterate) instead of switching to (or failing in) the feasibility problem
      [[nodiscard]] virtual bool terminate_at_line_search_failure(const Iterate& current_iterate, Iterate& trial_iterate,
         const Evaluations& current_evaluations, Evaluations& trial_evaluations) = 0;

      // threshold theta_min of the globalization strategy of the current phase
      [[nodiscard]] virtual double small_infeasibility_threshold() const = 0;

      [[nodiscard]] virtual std::string get_name() const = 0;
      [[nodiscard]] size_t get_number_subproblems_solved() const;
      // status of the last iterate tested for termination wrt the acceptable conditions (NOT_OPTIMAL if not acceptable)
      [[nodiscard]] SolutionStatus get_last_iterate_acceptable_status() const { return this->last_iterate_acceptable_status; }
      // number of updates of the parameterization (e.g. the barrier parameter) so far. Lets the globalization mechanism
      // detect that the problem definition changed (e.g. to reset the watchdog, as IPOPT does when mu changes)
      [[nodiscard]] size_t get_number_parameterization_updates() const { return this->number_parameterization_updates; }

   protected:
      const Norm progress_norm;
      const Norm residual_norm;
      const double primal_tolerance;
      const double dual_tolerance;
      const double diverging_iterate_threshold;
      const double unbounded_objective_threshold;
      size_t number_subproblems_solved{0};
      size_t number_parameterization_updates{0};

      // tight unscaled tolerances
      const double tight_unscaled_stationarity_tolerance;
      const double tight_unscaled_primal_infeasibility_tolerance;
      const double tight_unscaled_complementarity_tolerance;

      // acceptable termination (IPOPT-like)
      const double acceptable_tolerance;                    // scaled overall error [1e-6]
      const double acceptable_unscaled_stationarity_tolerance; // unscaled [1e10]
      const double acceptable_unscaled_primal_infeasibility_tolerance; // unscaled [1e-2]
      const double acceptable_unscaled_complementarity_tolerance;    // unscaled [1e-2]
      const double acceptable_objective_change_tolerance;   // relative [1e20]
      const size_t acceptable_iteration_threshold;          // [15], 0 disables the heuristic
      size_t acceptable_consecutive_iterations{0};
      std::optional<double> previous_objective{};
      // last iterate that satisfied the acceptable conditions (fallback upon line-search failure, see IPOPT)
      std::optional<Iterate> acceptable_iterate{};
      // last iterate of the current feasibility phase that was an acceptable infeasible stationary point
      std::optional<Iterate> acceptable_infeasible_iterate{};
      SolutionStatus last_iterate_acceptable_status{SolutionStatus::NOT_OPTIMAL};

      void compute_residuals(const OptimizationProblem& problem, Iterate& iterate, Evaluations& evaluations) const;
      [[nodiscard]] static Decision is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const OptimizationProblem& problem, const Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         Evaluations& trial_evaluations, const ProgressMeasures& predicted_reductions);
      static void accept_iterate_unconditionally(const OptimizationProblem& problem, Iterate& trial_iterate,
         Evaluations& trial_evaluations);
      [[nodiscard]] static bool satisfies_unscaled_tolerances(const OptimizationProblem& problem, const Iterate& iterate,
         double unscaled_primal_infeasibility, double stationarity_tolerance, double primal_infeasibility_tolerance,
         double complementarity_tolerance);
      [[nodiscard]] SolutionStatus check_termination(const OptimizationProblem& problem, const Iterate& iterate, Evaluations& evaluations);
      [[nodiscard]] SolutionStatus check_acceptable_termination(const OptimizationProblem& problem, const Iterate& iterate,
         const Evaluations& evaluations) const;
      [[nodiscard]] bool fall_back_after_line_search_failure(const OptimizationProblem& problem, const Iterate& current_iterate,
         Iterate& trial_iterate, const Evaluations& current_evaluations, Evaluations& trial_evaluations) const;
      [[nodiscard]] bool fall_back_to_infeasible_stationary_point(const OptimizationProblem& problem, const Iterate& current_iterate,
         Iterate& trial_iterate, const Evaluations& current_evaluations, Evaluations& trial_evaluations) const;
   };
} // namespace

#endif //UNO_CONSTRAINTRELAXATIONSTRATEGY_H