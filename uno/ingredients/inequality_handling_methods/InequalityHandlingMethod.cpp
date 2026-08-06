// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "InequalityHandlingMethod.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   bool InequalityHandlingMethod::is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const Subproblem& subproblem, const SolverWorkspace& solver_workspace, Iterate& current_iterate, Iterate& trial_iterate,
         const Direction& direction, double step_length, Evaluations& current_evaluations, Evaluations& trial_evaluations) const {
      this->problem.postprocess_iterate(trial_iterate);
      const double objective_multiplier = subproblem.problem.get_objective_multiplier();

      // evaluate progress measures
      this->evaluate_progress_measures(trial_iterate, trial_evaluations);
      trial_iterate.objective_multiplier = objective_multiplier;

      bool accept_iterate = false;
      if (direction.norm == 0.) {
         DEBUG << "Zero step acceptable\n";
         trial_evaluations.objective = this->problem.model.evaluate_objective(trial_iterate.primals);
         accept_iterate = true;
         statistics.set("Status", "0 primal step");
      }
      else {
         // determine acceptance wrt the globalization strategy
         Norm progress_norm = Norm::L1;
         const ProgressMeasures predicted_reductions = subproblem.compute_predicted_reductions(current_iterate, direction,
            step_length, progress_norm, current_evaluations, solver_workspace);
         accept_iterate = globalization_strategy.is_iterate_acceptable(statistics, current_iterate.progress, trial_iterate.progress,
            predicted_reductions, objective_multiplier);
         // check that the derivatives exist at the accepted trial iterate (an exception is thrown upon evaluation failure)
         if (accept_iterate) {
            trial_evaluations.evaluate_objective_gradient(this->problem.model, trial_iterate.primals);
            trial_evaluations.evaluate_jacobian(this->problem.model, trial_iterate.primals);
         }
      }
      return accept_iterate;
   }

   void InequalityHandlingMethod::evaluate_progress_measures(Iterate& iterate, Evaluations& evaluations) const {
      Norm progress_norm = Norm::L1;
      this->problem.set_infeasibility_measure(iterate, evaluations, progress_norm);
      this->problem.set_objective_measure(iterate, evaluations);
      this->problem.set_auxiliary_measure(iterate);
   }
} // namespace