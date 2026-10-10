// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_CONSTRAINTRELAXATIONSTRATEGY_H
#define UNO_CONSTRAINTRELAXATIONSTRATEGY_H

#include <cstddef>
#include <vector>
#include "ingredients/globalization_strategies/PredictedReductionModels.hpp"
#include "linear_algebra/Norm.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/SolutionStatus.hpp"
#include "options/Options.hpp"

namespace uno {
   // forward declarations
   class Direction;
   class EvaluationCache;
   class Evaluations;
   class GlobalizationStrategy;
   class Model;
   class Multipliers;
   class OptimizationProblem;
   class SolverWorkspace;
   class Statistics;
   class Subproblem;
   class UserCallbacks;
   class WarmstartInformation;

   class ConstraintRelaxationStrategy {
   public:
      explicit ConstraintRelaxationStrategy(const Options& options);
      virtual ~ConstraintRelaxationStrategy();

      virtual void initialize(Statistics& statistics, Iterate& initial_iterate, bool uses_trust_region,
         EvaluationCache& evaluation_cache, const Options& options, std::vector<OptionOverride>& option_overrides) = 0;

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
      virtual const Direction& compute_second_order_correction(Iterate& current_iterate) = 0;
      virtual void update_second_order_corrections(const Iterate& trial_iterate, Evaluations& trial_evaluations) = 0;

      // trial iterate acceptance
      [[nodiscard]] virtual PredictedReductionModels build_predicted_reduction_models(const Iterate& current_iterate,
         const Direction& direction, Evaluations& current_evaluations) const = 0;
      [[nodiscard]] virtual bool is_iterate_acceptable(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Iterate& trial_iterate, const Direction& direction, double step_length, bool uses_trust_region, Evaluations& current_evaluations,
         Evaluations& trial_evaluations, const ProgressMeasures& predicted_reductions, WarmstartInformation& warmstart_information,
         UserCallbacks& user_callbacks) = 0;
      // accept a tiny trial step without globalization test (IPOPT, Section 3.9)
      virtual void accept_tiny_step(Statistics& statistics, const Model& model, Iterate& current_iterate, Iterate& trial_iterate,
         const Direction& direction, double step_length, bool force_parameterization_update, Evaluations& current_evaluations,
         Evaluations& trial_evaluations, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) = 0;

      [[nodiscard]] virtual std::string get_name() const = 0;
      [[nodiscard]] size_t get_number_subproblems_solved() const;

   protected:
      const Norm residual_norm;
      const double primal_tolerance;
      const double dual_tolerance;
      const double loose_primal_tolerance;
      const double loose_dual_tolerance;
      size_t loose_tolerance_consecutive_iterations{0};
      const size_t loose_tolerance_iteration_threshold;
      const double diverging_iterate_threshold;
      const double unbounded_objective_threshold;
      size_t number_subproblems_solved{0};

      // acceptable termination (IPOPT-like)
      const double acceptable_tolerance;                    // scaled overall error [1e-6]
      const double acceptable_dual_infeasibility_tolerance; // unscaled [1e10]
      const double acceptable_primal_infeasibility_tolerance; // unscaled [1e-2]
      const double acceptable_complementarity_tolerance;    // unscaled [1e-2]
      const double acceptable_objective_change_tolerance;   // relative [1e20]
      const size_t acceptable_iteration_threshold;          // [15], 0 disables the heuristic
      size_t acceptable_consecutive_iterations{0};
      std::optional<double> previous_objective{};

      template <typename Problem>
      [[nodiscard]] SolutionStatus check_termination(const Problem& problem, Iterate& iterate, const Evaluations& evaluations);
      template <typename Problem>
      [[nodiscard]] SolutionStatus check_acceptable_termination(const Problem& problem, const Iterate& iterate,
         const Evaluations& evaluations) const;
   };

   template <typename Problem>
   SolutionStatus ConstraintRelaxationStrategy::check_termination(const Problem& problem, Iterate& iterate,
         const Evaluations& evaluations) {
      // unbounded iterate or objective
      if (norm_inf(view(iterate.primals, 0, problem.get_number_original_variables())) > this->diverging_iterate_threshold) {
         return SolutionStatus::DIVERGING_ITERATE;
      }
      if (iterate.scaled_residuals.primal_infeasibility <= this->primal_tolerance && evaluations.is_objective_computed &&
            evaluations.objective < this->unbounded_objective_threshold) {
         return SolutionStatus::UNBOUNDED_OBJECTIVE;
      }

      // tight termination: test convergence wrt the tight tolerance
      const SolutionStatus status_tight_tolerance = problem.check_first_order_convergence(iterate, this->primal_tolerance,
         this->dual_tolerance);
      if (status_tight_tolerance != SolutionStatus::NOT_OPTIMAL) {
         return status_tight_tolerance;
      }

      // acceptable termination: all acceptable conditions must hold for enough consecutive iterations
      const SolutionStatus status_acceptable = this->check_acceptable_termination(problem, iterate, evaluations);
      if (evaluations.is_objective_computed) {
         this->previous_objective = evaluations.objective;
      }
      if (status_acceptable == SolutionStatus::NOT_OPTIMAL) {
         this->acceptable_consecutive_iterations = 0;
         return SolutionStatus::NOT_OPTIMAL;
      }
      ++this->acceptable_consecutive_iterations;
      if (0 < this->acceptable_iteration_threshold &&
            this->acceptable_iteration_threshold <= this->acceptable_consecutive_iterations) {
         return status_acceptable;
      }
      return SolutionStatus::NOT_OPTIMAL;
   }

   // returns the convergence status if the iterate is acceptable, NOT_OPTIMAL otherwise.
   // Stateless (apart from the previous objective): can be reused as a fallback when the algorithm fails
   template <typename Problem>
   SolutionStatus ConstraintRelaxationStrategy::check_acceptable_termination(const Problem& problem,
         const Iterate& iterate, const Evaluations& evaluations) const {
      // 1. scaled overall optimality error
      const SolutionStatus status = problem.check_first_order_convergence(iterate, this->acceptable_tolerance,
         this->acceptable_tolerance);
      if (status == SolutionStatus::NOT_OPTIMAL) {
         return SolutionStatus::NOT_OPTIMAL;
      }
      // 2-4. unscaled residuals
      if (this->acceptable_dual_infeasibility_tolerance < iterate.scaled_residuals.stationarity ||
            this->acceptable_primal_infeasibility_tolerance < iterate.scaled_residuals.primal_infeasibility ||
            this->acceptable_complementarity_tolerance < iterate.scaled_residuals.complementarity) {
         return SolutionStatus::NOT_OPTIMAL;
      }
      // 5. relative change of the objective (unavailable at the first iteration: skip the test)
      if (evaluations.is_objective_computed && this->previous_objective.has_value()) {
         const double objective_change = std::abs(evaluations.objective - *this->previous_objective) /
            std::max(1., std::abs(evaluations.objective));
         if (this->acceptable_objective_change_tolerance < objective_change) {
            return SolutionStatus::NOT_OPTIMAL;
         }
      }
      return status;
   }
} // namespace

#endif //UNO_CONSTRAINTRELAXATIONSTRATEGY_H