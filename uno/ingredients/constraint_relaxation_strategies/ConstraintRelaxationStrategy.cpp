// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cassert>
#include <cmath>
#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "options/Options.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   ConstraintRelaxationStrategy::ConstraintRelaxationStrategy(const Options& options):
         progress_norm(norm_from_string(options.get_string("progress_norm"))),
         residual_norm(norm_from_string(options.get_string("residual_norm"))),
         primal_tolerance(options.get_double("primal_tolerance")),
         dual_tolerance(options.get_double("dual_tolerance")),
         diverging_iterate_threshold(options.get_double("diverging_iterate_threshold")),
         unbounded_objective_threshold(options.get_double("unbounded_objective_threshold")),
         tight_unscaled_stationarity_tolerance(options.get_double("tight_unscaled_stationarity_tolerance")),
         tight_unscaled_primal_infeasibility_tolerance(options.get_double("tight_unscaled_primal_infeasibility_tolerance")),
         tight_unscaled_complementarity_tolerance(options.get_double("tight_unscaled_complementarity_tolerance")),
         acceptable_tolerance(options.get_double("acceptable_tolerance")),
         acceptable_unscaled_stationarity_tolerance(options.get_double("acceptable_unscaled_stationarity_tolerance")),
         acceptable_unscaled_primal_infeasibility_tolerance(options.get_double("acceptable_unscaled_primal_infeasibility_tolerance")),
         acceptable_unscaled_complementarity_tolerance(options.get_double("acceptable_unscaled_complementarity_tolerance")),
         acceptable_objective_change_tolerance(options.get_double("acceptable_objective_change_tolerance")),
         acceptable_iteration_threshold(options.get_unsigned_int("acceptable_iteration_threshold")) {
   }

   ConstraintRelaxationStrategy::~ConstraintRelaxationStrategy() = default;

   // stationarity errors:
   // - for KKT conditions: with standard multipliers and current objective multiplier
   // - for FJ conditions: with standard multipliers and 0 objective multiplier
   // - for feasibility problem: with feasibility multipliers and 0 objective multiplier
   void ConstraintRelaxationStrategy::compute_residuals(const OptimizationProblem& problem, Iterate& iterate,
         Evaluations& evaluations) const {
      // stationarity error (norm of the Lagrangian gradient)
      problem.evaluate_lagrangian_gradient(iterate.primals, iterate.multipliers.constraints, iterate.multipliers.lower_bounds,
         iterate.multipliers.upper_bounds, evaluations, iterate.scaled_residuals.lagrangian_gradient.view());
      iterate.scaled_residuals.stationarity = norm(this->residual_norm, iterate.scaled_residuals.lagrangian_gradient);

      // primal infeasibility (constraint violation of the model, possibly with additional residuals of the reformulation)
      iterate.scaled_residuals.primal_infeasibility = problem.primal_infeasibility(iterate, evaluations, this->residual_norm);
      assert(evaluations.are_constraints_computed);

      // complementarity error
      iterate.scaled_residuals.complementarity = problem.complementarity_error(iterate.primals, evaluations.constraints,
         iterate.multipliers, this->residual_norm);

      // scaling factors
      iterate.scaled_residuals.stationarity_scaling = problem.compute_stationarity_scaling(iterate.multipliers);
      iterate.scaled_residuals.complementarity_scaling = problem.compute_complementarity_scaling(iterate.multipliers);
   }

   Decision ConstraintRelaxationStrategy::is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const OptimizationProblem& problem, const Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         Evaluations& trial_evaluations, const ProgressMeasures& predicted_reductions) {
      problem.postprocess_iterate(trial_iterate, trial_evaluations);
      const double objective_multiplier = problem.get_objective_multiplier();
      trial_iterate.objective_multiplier = objective_multiplier;

      // zero step in all the primal variables (direction.norm only measures the original variables: slacks or elastic
      // variables may still move)
      if (norm_inf(direction.primals) == 0.) {
         DEBUG << "Zero step acceptable\n";
         trial_evaluations.evaluate_objective(problem.model, trial_iterate.primals);
         statistics.set("Status", "0 primal step");
         return Decision::accept();
      }
      if (is_infinite(trial_iterate.progress.auxiliary)) {
         DEBUG << "The auxiliary measure is infinite (iterate too close to the bounds), rejecting the trial iterate\n";
         statistics.set("Status", "inf auxiliary");
         return Decision::reject(RejectionReason::INVALID_TRIAL_ITERATE);
      }
      // determine acceptance wrt the globalization strategy
      const Decision decision = globalization_strategy.is_iterate_acceptable(statistics, current_iterate.progress,
         trial_iterate.progress, predicted_reductions, objective_multiplier);
      if (decision.accepted) {
         // check that the derivatives exist at the accepted trial iterate (an exception is thrown upon evaluation failure)
         trial_evaluations.evaluate_objective_gradient(problem.model, trial_iterate.primals);
         trial_evaluations.evaluate_jacobian(problem.model, trial_iterate.primals);
      }
      return decision;
   }

   void ConstraintRelaxationStrategy::accept_iterate_unconditionally(const OptimizationProblem& problem, Iterate& trial_iterate,
         Evaluations& trial_evaluations) {
      problem.postprocess_iterate(trial_iterate, trial_evaluations);
      trial_iterate.objective_multiplier = problem.get_objective_multiplier();
      // the next subproblem must be well defined (an exception is thrown upon evaluation failure)
      trial_evaluations.evaluate_objective_gradient(problem.model, trial_iterate.primals);
      trial_evaluations.evaluate_jacobian(problem.model, trial_iterate.primals);
   }

   bool ConstraintRelaxationStrategy::satisfies_unscaled_tolerances(const OptimizationProblem& problem, const Iterate& iterate,
         double unscaled_primal_infeasibility, double stationarity_tolerance, double primal_infeasibility_tolerance,
         double complementarity_tolerance) {
      const Model& model = problem.model;
      return model.unscaled_stationarity(iterate.scaled_residuals.stationarity) <= stationarity_tolerance &&
         unscaled_primal_infeasibility <= primal_infeasibility_tolerance &&
         model.unscaled_complementarity(iterate.scaled_residuals.complementarity) <= complementarity_tolerance;
   }

   SolutionStatus ConstraintRelaxationStrategy::check_termination(const OptimizationProblem& problem, const Iterate& iterate,
         Evaluations& evaluations) {
      // diverging iterate
      if (norm_inf(view(iterate.primals, 0, problem.get_number_original_variables())) > this->diverging_iterate_threshold) {
         return SolutionStatus::DIVERGING_ITERATE;
      }
      // unbounded objective + feasible
      evaluations.evaluate_constraints(problem.model, iterate.primals); // no-op if cached or unconstrained
      const double unscaled_primal_infeasibility = problem.model.unscaled_primal_infeasibility(iterate.scaled_residuals.primal_infeasibility,
         evaluations.constraints, this->residual_norm);
      if (unscaled_primal_infeasibility <= this->tight_unscaled_primal_infeasibility_tolerance) {
         evaluations.evaluate_objective(problem.model, iterate.primals);
         if (problem.model.unscaled_objective(evaluations.objective) < this->unbounded_objective_threshold) {
            return SolutionStatus::UNBOUNDED_OBJECTIVE;
         }
      }

      // tight termination: test convergence wrt the tight tolerance
      const SolutionStatus status_tight_tolerance = problem.check_first_order_convergence(iterate, this->primal_tolerance,
         this->dual_tolerance);
      if (status_tight_tolerance == SolutionStatus::FEASIBLE_KKT_POINT && satisfies_unscaled_tolerances(problem, iterate,
            unscaled_primal_infeasibility, this->tight_unscaled_stationarity_tolerance,
            this->tight_unscaled_primal_infeasibility_tolerance, this->tight_unscaled_complementarity_tolerance)) {
         return status_tight_tolerance;
      }
      else if (status_tight_tolerance == SolutionStatus::INFEASIBLE_STATIONARY_POINT) {
         // the point should also be infeasible for the original (unscaled) constraints
         if (this->tight_unscaled_primal_infeasibility_tolerance < unscaled_primal_infeasibility) {
            return status_tight_tolerance;
         }
      }

      // acceptable termination: all acceptable conditions must hold for enough consecutive iterations
      const SolutionStatus status_acceptable = (0 < this->acceptable_iteration_threshold) ?
         this->check_acceptable_termination(problem, iterate, evaluations) : SolutionStatus::NOT_OPTIMAL;
      if (evaluations.is_objective_computed) {
         this->previous_objective = evaluations.objective; // after the objective-change test
      }
      this->last_iterate_acceptable_status = status_acceptable;
      if (status_acceptable == SolutionStatus::NOT_OPTIMAL) {
         this->acceptable_consecutive_iterations = 0;
         return SolutionStatus::NOT_OPTIMAL;
      }
      if (status_acceptable == SolutionStatus::INFEASIBLE_STATIONARY_POINT) {
         this->acceptable_infeasible_iterate = iterate;
      }
      else {
         this->acceptable_iterate = iterate;
      }
      ++this->acceptable_consecutive_iterations;
      if (this->acceptable_iteration_threshold <= this->acceptable_consecutive_iterations) {
         return status_acceptable;
      }
      return SolutionStatus::NOT_OPTIMAL;
   }

   // returns the convergence status if the iterate is acceptable, NOT_OPTIMAL otherwise.
   // Stateless (apart from the previous objective): can be reused as a fallback when the algorithm fails
   SolutionStatus ConstraintRelaxationStrategy::check_acceptable_termination(const OptimizationProblem& problem,
         const Iterate& iterate, const Evaluations& evaluations) const {
      // 1. scaled overall optimality error
      const SolutionStatus status = problem.check_first_order_convergence(iterate, this->acceptable_tolerance,
         this->acceptable_tolerance);
      const double unscaled_primal_infeasibility = problem.model.unscaled_primal_infeasibility(iterate.scaled_residuals.primal_infeasibility,
         evaluations.constraints, this->residual_norm);
      // acceptable infeasible stationary point (IPOPT: the restoration phase converged to an acceptable point of local
      // infeasibility). As in the tight test, the point should also be infeasible for the original (unscaled) constraints
      if (status == SolutionStatus::INFEASIBLE_STATIONARY_POINT) {
         return (this->tight_unscaled_primal_infeasibility_tolerance < unscaled_primal_infeasibility) ? status :
            SolutionStatus::NOT_OPTIMAL;
      }
      if (status != SolutionStatus::FEASIBLE_KKT_POINT) {
         return SolutionStatus::NOT_OPTIMAL;
      }
      // 2-4. unscaled residuals
      if (!satisfies_unscaled_tolerances(problem, iterate, unscaled_primal_infeasibility, this->acceptable_unscaled_stationarity_tolerance,
            this->acceptable_unscaled_primal_infeasibility_tolerance, this->acceptable_unscaled_complementarity_tolerance)) {
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

   // IPOPT (BacktrackingLineSearch::FindAcceptableTrialPoint), upon line-search failure in the optimality phase:
   // - if the current iterate is acceptable, terminate there;
   // - if the current iterate is almost feasible, the feasibility restoration cannot make progress: terminate at the last
   //   acceptable iterate if any, otherwise fail.
   // Returns true if the trial iterate was set to a terminal iterate
   bool ConstraintRelaxationStrategy::fall_back_after_line_search_failure(const OptimizationProblem& problem, const Iterate& current_iterate,
         Iterate& trial_iterate, const Evaluations& current_evaluations, Evaluations& trial_evaluations) const {
      if (0 < this->acceptable_iteration_threshold && this->check_acceptable_termination(problem, current_iterate,
            current_evaluations) == SolutionStatus::FEASIBLE_KKT_POINT) {
         WARNING << "The line search failed at an acceptable iterate\n";
         trial_iterate = current_iterate;
         trial_evaluations = current_evaluations;
         trial_iterate.status = SolutionStatus::FEASIBLE_KKT_POINT;
         return true;
      }
      const double unscaled_primal_infeasibility = problem.model.unscaled_primal_infeasibility(
         current_iterate.scaled_residuals.primal_infeasibility, current_evaluations.constraints, Norm::INF);
      if (current_iterate.progress.infeasibility <= 1e-2 * this->primal_tolerance &&
            unscaled_primal_infeasibility <= 1e-1 * this->tight_unscaled_primal_infeasibility_tolerance) {
         if (this->acceptable_iterate.has_value()) {
            WARNING << "The line search failed at an almost feasible point, the last acceptable iterate is restored\n";
            trial_iterate = *this->acceptable_iterate;
            trial_evaluations.reset(); // the functions are evaluated again at the restored iterate upon postprocessing
            trial_iterate.status = SolutionStatus::FEASIBLE_KKT_POINT;
            return true;
         }
         throw std::runtime_error("The line search failed at an almost feasible point and no acceptable iterate is available");
      }
      return false;
   }

   // upon line-search failure in the feasibility phase: terminate at the current iterate if it is an acceptable infeasible
   // stationary point, otherwise at the last such iterate of the current feasibility phase (if any).
   // Returns true if the trial iterate was set to a terminal iterate
   bool ConstraintRelaxationStrategy::fall_back_to_infeasible_stationary_point(const OptimizationProblem& problem,
         const Iterate& current_iterate, Iterate& trial_iterate, const Evaluations& current_evaluations,
         Evaluations& trial_evaluations) const {
      if (this->acceptable_iteration_threshold == 0) {
         return false;
      }
      if (this->check_acceptable_termination(problem, current_iterate, current_evaluations) ==
            SolutionStatus::INFEASIBLE_STATIONARY_POINT) {
         WARNING << "The line search failed at an acceptable infeasible stationary point\n";
         trial_iterate = current_iterate;
         trial_evaluations = current_evaluations;
            }
      else if (this->acceptable_infeasible_iterate.has_value()) {
         WARNING << "The line search failed in the feasibility phase, the last acceptable infeasible stationary point is restored\n";
         trial_iterate = *this->acceptable_infeasible_iterate;
         trial_evaluations.reset(); // the functions are evaluated again at the restored iterate upon postprocessing
      }
      else {
         return false;
      }
      trial_iterate.status = SolutionStatus::INFEASIBLE_STATIONARY_POINT;
      return true;
   }

   size_t ConstraintRelaxationStrategy::get_number_subproblems_solved() const {
      return this->number_subproblems_solved;
   }
} // namespace