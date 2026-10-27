// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "InequalityHandlingMethod.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   InequalityHandlingMethod::InequalityHandlingMethod(const OptimizationProblem& problem, const Options& options):
         problem(problem),
         progress_norm(norm_from_string(options.get_string("progress_norm"))),
         residual_norm(norm_from_string(options.get_string("residual_norm"))) {
   }

   // protected member functions

   void InequalityHandlingMethod::evaluate_progress_measures(const OptimizationProblem& problem, Iterate& iterate,
         Evaluations& evaluations) const {
      problem.set_infeasibility_measure(iterate, evaluations, this->progress_norm);
      problem.set_objective_measure(iterate, evaluations);
      problem.set_auxiliary_measure(iterate);
   }

   // stationarity errors:
   // - for KKT conditions: with standard multipliers and current objective multiplier
   // - for FJ conditions: with standard multipliers and 0 objective multiplier
   // - for feasibility problem: with feasibility multipliers and 0 objective multiplier
   void InequalityHandlingMethod::compute_residuals(const OptimizationProblem& problem, Iterate& iterate,
         Evaluations& evaluations) const {
      // stationarity error (norm of the Lagrangian gradient)
      problem.evaluate_lagrangian_gradient(iterate, evaluations, iterate.scaled_residuals.lagrangian_gradient);
      iterate.scaled_residuals.stationarity = norm(this->residual_norm, iterate.scaled_residuals.lagrangian_gradient);

      // primal infeasibility (constraint violation of the model, possibly with additional residuals of the reformulation)
      iterate.scaled_residuals.primal_infeasibility = problem.primal_infeasibility(iterate, evaluations, this->residual_norm);
      evaluations.evaluate_constraints(problem.model, iterate.primals); // no-op if cached

      // complementarity error
      iterate.scaled_residuals.complementarity = problem.complementarity_error(iterate.primals, evaluations.constraints,
         iterate.multipliers, this->residual_norm);

      // scaling factors
      iterate.scaled_residuals.stationarity_scaling = problem.compute_stationarity_scaling(iterate.multipliers);
      iterate.scaled_residuals.complementarity_scaling = problem.compute_complementarity_scaling(iterate.multipliers);
   }

   Decision InequalityHandlingMethod::is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const Subproblem& subproblem, const Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         Evaluations& trial_evaluations, const ProgressMeasures& predicted_reductions) const {
      subproblem.problem.postprocess_iterate(trial_iterate, trial_evaluations);
      const double objective_multiplier = subproblem.problem.get_objective_multiplier();

      // evaluate progress measures
      evaluate_progress_measures(subproblem.problem, trial_iterate, trial_evaluations);
      trial_iterate.objective_multiplier = objective_multiplier;

      // zero step in all the primal variables (direction.norm only measures the original variables: slacks or elastic
      // variables may still move)
      if (norm_inf(direction.primals) == 0.) {
         DEBUG << "Zero step acceptable\n";
         trial_evaluations.evaluate_objective(this->problem.model, trial_iterate.primals);
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
         trial_evaluations.evaluate_objective_gradient(this->problem.model, trial_iterate.primals);
         trial_evaluations.evaluate_jacobian(this->problem.model, trial_iterate.primals);
      }
      return decision;
   }

   void InequalityHandlingMethod::accept_iterate_unconditionally(const Subproblem& subproblem, Iterate& trial_iterate,
         Evaluations& trial_evaluations) const {
      subproblem.problem.postprocess_iterate(trial_iterate, trial_evaluations);
      evaluate_progress_measures(subproblem.problem, trial_iterate, trial_evaluations);
      trial_iterate.objective_multiplier = subproblem.problem.get_objective_multiplier();
      // the next subproblem must be well defined (an exception is thrown upon evaluation failure)
      trial_evaluations.evaluate_objective_gradient(this->problem.model, trial_iterate.primals);
      trial_evaluations.evaluate_jacobian(this->problem.model, trial_iterate.primals);
   }
} // namespace