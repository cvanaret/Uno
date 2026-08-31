// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_CONSTRAINTRELAXATIONSTRATEGY_H
#define UNO_CONSTRAINTRELAXATIONSTRATEGY_H

#include <cstddef>
#include "linear_algebra/Norm.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/SolutionStatus.hpp"

namespace uno {
   // forward declarations
   class Direction;
   class EvaluationCache;
   class Evaluations;
   class GlobalizationStrategy;
   class Model;
   class Multipliers;
   class OptimizationProblem;
   class Options;
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
         EvaluationCache& evaluation_cache, Options& options) = 0;

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
      [[nodiscard]] virtual bool is_iterate_acceptable(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Iterate& trial_iterate, const Direction& direction, double step_length, bool uses_trust_region,
         Evaluations& current_evaluations, Evaluations& trial_evaluations, WarmstartInformation& warmstart_information,
         UserCallbacks& user_callbacks) = 0;

      [[nodiscard]] virtual std::string get_name() const = 0;
      [[nodiscard]] size_t get_number_subproblems_solved() const;

   protected:
      const Norm residual_norm;
      const double residual_scaling_threshold;
      const double primal_tolerance;
      const double dual_tolerance;
      const double loose_primal_tolerance;
      const double loose_dual_tolerance;
      size_t loose_tolerance_consecutive_iterations{0};
      const size_t loose_tolerance_iteration_threshold;
      const double diverging_iterate_threshold;
      const double unbounded_objective_threshold;
      size_t number_subproblems_solved{0};

      void compute_residuals(const OptimizationProblem& problem, Iterate& iterate, Evaluations& evaluations) const;
      [[nodiscard]] double compute_stationarity_scaling(const Model& model, const Multipliers& multipliers) const;
      [[nodiscard]] double compute_complementarity_scaling(const Model& model, const Multipliers& multipliers) const;

      template <typename Problem>
      [[nodiscard]] SolutionStatus check_termination(const Problem& problem, Iterate& iterate, const Evaluations& evaluations);
   };

   template <typename Problem>
   SolutionStatus ConstraintRelaxationStrategy::check_termination(const Problem& problem, Iterate& iterate,
         const Evaluations& evaluations) {
      // unbounded iterate or objective
      if (norm_inf(view(iterate.primals, 0, problem.get_number_original_variables())) > this->diverging_iterate_threshold) {
         return SolutionStatus::DIVERGING_ITERATE;
      }
      if (iterate.primal_infeasibility <= this->primal_tolerance && evaluations.is_objective_computed &&
            evaluations.objective < this->unbounded_objective_threshold) {
         return SolutionStatus::UNBOUNDED_OBJECTIVE;
      }

      // test convergence wrt the tight tolerance
      const SolutionStatus status_tight_tolerance = problem.check_first_order_convergence(iterate, this->primal_tolerance,
         this->dual_tolerance);
      if (status_tight_tolerance != SolutionStatus::NOT_OPTIMAL || (this->loose_primal_tolerance <= this->primal_tolerance &&
            this->loose_dual_tolerance <= this->dual_tolerance)) {
         return status_tight_tolerance;
      }

      // if not converged, check convergence wrt loose tolerances (provided they are strictly looser than the tight tolerances)
      const SolutionStatus status_loose_tolerance = problem.check_first_order_convergence(iterate, this->loose_primal_tolerance,
         this->loose_dual_tolerance);
      // if converged, keep track of the number of consecutive iterations
      if (status_loose_tolerance != SolutionStatus::NOT_OPTIMAL) {
         ++this->loose_tolerance_consecutive_iterations;
      }
      else {
         return SolutionStatus::NOT_OPTIMAL;
      }
      // check if loose tolerance achieved for enough consecutive iterations
      if (this->loose_tolerance_iteration_threshold <= this->loose_tolerance_consecutive_iterations) {
         return status_loose_tolerance;
      }
      else {
         return SolutionStatus::NOT_OPTIMAL;
      }
   }
} // namespace

#endif //UNO_CONSTRAINTRELAXATIONSTRATEGY_H