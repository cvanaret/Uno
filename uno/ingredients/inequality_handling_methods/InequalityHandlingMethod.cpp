// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "InequalityHandlingMethod.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "options/Options.hpp"

namespace uno {
   InequalityHandlingMethod::InequalityHandlingMethod(const OptimizationProblem& problem, const Options& options):
         problem(problem),
         progress_norm(norm_from_string(options.get_string("progress_norm"))),
         residual_norm(norm_from_string(options.get_string("residual_norm"))) {
   }

   // protected member functions

   void InequalityHandlingMethod::accept_iterate_unconditionally(const Subproblem& subproblem, Iterate& trial_iterate,
         Evaluations& trial_evaluations) const {
      subproblem.problem.postprocess_iterate(trial_iterate, trial_evaluations);
      trial_iterate.objective_multiplier = subproblem.problem.get_objective_multiplier();
      // the next subproblem must be well defined (an exception is thrown upon evaluation failure)
      trial_evaluations.evaluate_objective_gradient(this->problem.model, trial_iterate.primals);
      trial_evaluations.evaluate_jacobian(this->problem.model, trial_iterate.primals);
   }
} // namespace