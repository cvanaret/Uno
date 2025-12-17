// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "InequalityHandlingMethod.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"
#include "tools/UserCallbacks.hpp"

namespace uno {
   InequalityHandlingMethod::InequalityHandlingMethod(const Options& options):
      progress_norm(norm_from_string(options.get_string("progress_norm"))) {
   }

   void InequalityHandlingMethod::evaluate_progress_measures(const OptimizationProblem& problem, Iterate& iterate) const {
      problem.set_infeasibility_measure(iterate, this->progress_norm);
      problem.set_objective_measure(iterate);
      problem.set_auxiliary_measure(iterate);
   }

   bool InequalityHandlingMethod::is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const Subproblem& subproblem, const EvaluationSpace& evaluation_space, Iterate& current_iterate, Iterate& trial_iterate,
         const Direction& direction, double step_length, UserCallbacks& user_callbacks) {
      this->postprocess_iterate(trial_iterate);
      const double objective_multiplier = subproblem.problem.get_objective_multiplier();

      // evaluate progress measures
      trial_iterate.objective_multiplier = objective_multiplier;
      if (this->subproblem_definition_changed) {
         DEBUG << "The subproblem definition changed, the globalization strategy is reset and the auxiliary measure is recomputed\n";
         globalization_strategy.reset();
         subproblem.problem.set_auxiliary_measure(current_iterate);
         this->subproblem_definition_changed = false;
      }
      this->evaluate_progress_measures(subproblem.problem, trial_iterate);

      bool accept_iterate = false;
      if (direction.norm == 0.) {
         DEBUG << "Zero step acceptable\n";
         trial_iterate.evaluate_objective(subproblem.problem.model);
         accept_iterate = true;
         statistics.set("Status", "0 primal step");
      }
      else {
         const ProgressMeasures predicted_reductions = subproblem.compute_predicted_reductions(evaluation_space, direction,
            step_length, this->progress_norm);
         accept_iterate = globalization_strategy.is_iterate_acceptable(statistics, current_iterate.progress, trial_iterate.progress,
            predicted_reductions, objective_multiplier);
      }
      if (accept_iterate) {
         user_callbacks.notify_acceptable_iterate(trial_iterate.primals, trial_iterate.multipliers, objective_multiplier,
            trial_iterate.progress.infeasibility, trial_iterate.residuals.stationarity, trial_iterate.residuals.complementarity);
      }
      return accept_iterate;
   }
} // namespace