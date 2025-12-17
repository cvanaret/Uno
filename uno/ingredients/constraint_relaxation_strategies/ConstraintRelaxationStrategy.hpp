// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_CONSTRAINTRELAXATIONSTRATEGY_H
#define UNO_CONSTRAINTRELAXATIONSTRATEGY_H

#include <cstddef>
#include "linear_algebra/Norm.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/SolutionStatus.hpp"

namespace uno {
   // forward declarations
   class Direction;
   class Model;
   class Multipliers;
   class OptimizationProblem;
   class Options;
   class Statistics;
   class UserCallbacks;
   class WarmstartInformation;

   class ConstraintRelaxationStrategy {
   public:
      explicit ConstraintRelaxationStrategy(const Options& options);
      virtual ~ConstraintRelaxationStrategy();

      virtual void initialize(Statistics& statistics, const Model& model, Iterate& initial_iterate, Direction& direction,
         double trust_region_radius) = 0;

      // direction computation
      virtual void compute_feasible_direction(Statistics& statistics, Iterate& current_iterate, Direction& direction,
         double trust_region_radius, WarmstartInformation& warmstart_information) = 0;
      [[nodiscard]] virtual bool solving_feasibility_problem() const = 0;
      virtual void switch_to_feasibility_problem(Statistics& statistics, Iterate& current_iterate, double trust_region_radius,
         WarmstartInformation& warmstart_information) = 0;

      // trial iterate acceptance
      [[nodiscard]] virtual bool is_iterate_acceptable(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Iterate& trial_iterate, const Direction& direction,
         double step_length, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) = 0;
      [[nodiscard]] virtual SolutionStatus check_termination(const Model& model, Iterate& iterate) = 0;

      [[nodiscard]] virtual std::string get_name() const = 0;
      [[nodiscard]] virtual size_t get_number_subproblems_solved() const = 0;

   protected:
      const Norm residual_norm;
      const double residual_scaling_threshold;
      const double primal_tolerance;
      const double dual_tolerance;
      const double loose_primal_tolerance;
      const double loose_dual_tolerance;
      size_t loose_tolerance_consecutive_iterations{0};
      const size_t loose_tolerance_consecutive_iteration_threshold;
      const double unbounded_objective_threshold;

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
      if (status_tight_tolerance != SolutionStatus::NOT_OPTIMAL || (this->loose_primal_tolerance <= this->primal_tolerance &&
            this->loose_dual_tolerance <= this->dual_tolerance)) {
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