// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_CONSTRAINTRELAXATIONSTRATEGY_H
#define UNO_CONSTRAINTRELAXATIONSTRATEGY_H

#include <cstddef>
#include <functional>
#include "ingredients/globalization_strategies/ProgressMeasures.hpp"
#include "linear_algebra/Norm.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/SolutionStatus.hpp"

namespace uno {
   // forward declarations
   class Direction;
   class GlobalizationStrategy;
   class InequalityHandlingMethod;
   class Model;
   class Multipliers;
   class OptimizationProblem;
   class Options;
   class Statistics;
   class UserCallbacks;
   template <typename ElementType>
   class Vector;
   class WarmstartInformation;

   class ConstraintRelaxationStrategy {
   public:
      explicit ConstraintRelaxationStrategy(const Options& options);
      virtual ~ConstraintRelaxationStrategy();

      virtual void initialize(Statistics& statistics, const Model& model, Iterate& initial_iterate, Direction& direction,
         double trust_region_radius, const Options& options) = 0;

      // direction computation
      virtual void compute_feasible_direction(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const Model& model, Iterate& current_iterate, Direction& direction, double trust_region_radius,
         WarmstartInformation& warmstart_information) = 0;
      [[nodiscard]] virtual bool solving_feasibility_problem() const = 0;
      virtual void switch_to_feasibility_problem(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const Model& model, Iterate& current_iterate, double trust_region_radius, WarmstartInformation& warmstart_information) = 0;

      // trial iterate acceptance
      [[nodiscard]] virtual bool is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const Model& model, Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, double step_length,
         WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) = 0;
      [[nodiscard]] virtual SolutionStatus check_termination(const Model& model, Iterate& iterate) = 0;

      [[nodiscard]] virtual std::string get_name() const = 0;
      [[nodiscard]] virtual size_t get_hessian_evaluation_count() const = 0;
      [[nodiscard]] virtual size_t get_number_subproblems_solved() const = 0;

   protected:
      const Norm progress_norm;
      const Norm residual_norm;
      const double residual_scaling_threshold;
      const double primal_tolerance;
      const double dual_tolerance;
      const double loose_dual_tolerance;
      size_t loose_tolerance_consecutive_iterations{0};
      const size_t loose_tolerance_consecutive_iteration_threshold;
      const double unbounded_objective_threshold;
      // first_order_predicted_reduction is true when the predicted reduction can be taken as first-order (e.g. in line-search methods)
      const bool first_order_predicted_reduction;

      void set_objective_measure(const Model& model, Iterate& iterate) const;
      void set_infeasibility_measure(const Model& model, Iterate& iterate) const;
      [[nodiscard]] double compute_predicted_infeasibility_reduction(InequalityHandlingMethod& inequality_handling_method,
         const Model& model, const Iterate& current_iterate, const Vector<double>& primal_direction, double step_length) const;
      [[nodiscard]] std::function<double(double)> compute_predicted_objective_reduction(InequalityHandlingMethod& inequality_handling_method,
         const Iterate& current_iterate, const Vector<double>& primal_direction, double step_length) const;
      void compute_progress_measures(InequalityHandlingMethod& inequality_handling_method, const OptimizationProblem& problem,
         GlobalizationStrategy& globalization_strategy, Iterate& current_iterate, Iterate& trial_iterate) const;
      [[nodiscard]] ProgressMeasures compute_predicted_reductions(InequalityHandlingMethod& inequality_handling_method,
         const OptimizationProblem& problem, const Iterate& current_iterate, const Direction& direction, double step_length) const;
      [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const OptimizationProblem& problem, InequalityHandlingMethod& inequality_handling_method, Iterate& current_iterate,
         Iterate& trial_iterate, const Direction& direction, double step_length, UserCallbacks& user_callbacks) const;
      virtual void evaluate_progress_measures(InequalityHandlingMethod& inequality_handling_method,
         const OptimizationProblem& problem, Iterate& iterate) const = 0;

      void compute_primal_dual_residuals(const OptimizationProblem& problem, Iterate& iterate) const;
      [[nodiscard]] double compute_stationarity_scaling(const Model& model, const Multipliers& multipliers) const;
      [[nodiscard]] double compute_complementarity_scaling(const Model& model, const Multipliers& multipliers) const;

      template <typename Problem>
      [[nodiscard]] SolutionStatus check_termination(const Problem& problem, Iterate& iterate);
   };

   template <typename Problem>
   SolutionStatus ConstraintRelaxationStrategy::check_termination(const Problem& problem, Iterate& iterate) {
      if (iterate.is_objective_computed && iterate.evaluations.objective < this->unbounded_objective_threshold) {
         return SolutionStatus::UNBOUNDED;
      }

      // test convergence wrt the tight tolerance
      const SolutionStatus status_tight_tolerance = problem.check_first_order_convergence(iterate, this->primal_tolerance,
         this->dual_tolerance);
      if (status_tight_tolerance != SolutionStatus::NOT_OPTIMAL || this->loose_dual_tolerance <= this->primal_tolerance) {
         return status_tight_tolerance;
      }

      // if not converged, check convergence wrt loose tolerance (provided it is strictly looser than the tight tolerance)
      const SolutionStatus status_loose_tolerance = problem.check_first_order_convergence(iterate, this->primal_tolerance,
         this->loose_dual_tolerance);
      // if converged, keep track of the number of consecutive iterations
      if (status_loose_tolerance != SolutionStatus::NOT_OPTIMAL) {
         ++this->loose_tolerance_consecutive_iterations;
      }
      else {
         this->loose_tolerance_consecutive_iterations = 0;
         return SolutionStatus::NOT_OPTIMAL;
      }
      // check if loose tolerance achieved for enough consecutive iterations
      if (this->loose_tolerance_consecutive_iteration_threshold <= this->loose_tolerance_consecutive_iterations) {
         return status_loose_tolerance;
      }
      else {
         return SolutionStatus::NOT_OPTIMAL;
      }
   }
} // namespace

#endif //UNO_CONSTRAINTRELAXATIONSTRATEGY_H