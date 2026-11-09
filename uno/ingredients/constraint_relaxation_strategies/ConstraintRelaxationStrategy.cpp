// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ConstraintRelaxationStrategy.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "optimization/Direction.hpp"
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

   size_t ConstraintRelaxationStrategy::get_number_subproblems_solved() const {
      return this->number_subproblems_solved;
   }
} // namespace