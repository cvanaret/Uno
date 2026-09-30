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
#include "tools/Symbols.hpp"

namespace uno {
   InequalityHandlingMethod::InequalityHandlingMethod(const OptimizationProblem& problem, const Options& options):
         problem(problem), progress_norm(norm_from_string(options.get_string("progress_norm"))),
         theta_min(options.get_double("barrier_small_infeasibility_factor")) {
   }

   // protected member functions

   void InequalityHandlingMethod::evaluate_progress_measures(const OptimizationProblem& problem, Iterate& iterate,
         Evaluations& evaluations) const {
      problem.set_infeasibility_measure(iterate, evaluations, this->progress_norm);
      problem.set_objective_measure(iterate, evaluations);
      problem.set_auxiliary_measure(iterate);
   }

   bool InequalityHandlingMethod::is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         const Subproblem& subproblem, const Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction,
         Evaluations& trial_evaluations, const ProgressMeasures& predicted_reductions, bool is_full_step) const {
      subproblem.problem.postprocess_iterate(trial_iterate);
      const double objective_multiplier = subproblem.problem.get_objective_multiplier();

      // evaluate progress measures
      evaluate_progress_measures(subproblem.problem, trial_iterate, trial_evaluations);
      trial_iterate.objective_multiplier = objective_multiplier;

      const bool tiny_direction = is_tiny_direction(current_iterate, direction);
      if (!tiny_direction) {
         this->number_consecutive_tiny_directions = 0;
      }

      bool accept_iterate = false;
      if (direction.norm == 0.) {
         DEBUG << "Zero step acceptable\n";
         trial_evaluations.evaluate_objective(this->problem.model, trial_iterate.primals);
         accept_iterate = true;
         statistics.set("Status", "0 primal step");
      }
      else if (is_infinite(trial_iterate.progress.auxiliary)) {
         DEBUG << "The auxiliary measure is infinite (iterate too close to the bounds), rejecting the trial iterate\n";
         accept_iterate = false;
         statistics.set("Status", "inf auxiliary");
      }
      else {
         // determine acceptance wrt the globalization strategy
         accept_iterate = globalization_strategy.is_iterate_acceptable(statistics, current_iterate.progress, trial_iterate.progress,
            predicted_reductions, objective_multiplier);

         // if rejected, accept it if (full) tiny direction
         if (!accept_iterate && is_full_step && tiny_direction) {
            ++this->number_consecutive_tiny_directions;
            if (this->number_consecutive_tiny_directions >= this->consecutive_tiny_directions_threshold) {
               accept_iterate = true;
               DEBUG << "Accepting tiny step\n";
               statistics.set("Status", std::string(symbols::check) + " (tiny)");
               this->number_consecutive_tiny_directions = 0;
            }
         }
      }

      // last chance to reject the trial iterate: check that the functions and their derivatives exist
      // (an exception is thrown upon evaluation failure)
      if (accept_iterate) {
         trial_evaluations.evaluate_constraints(this->problem.model, trial_iterate.primals);
         trial_evaluations.evaluate_objective_gradient(this->problem.model, trial_iterate.primals);
         trial_evaluations.evaluate_jacobian(this->problem.model, trial_iterate.primals);
      }
      return accept_iterate;
   }

   bool InequalityHandlingMethod::is_tiny_direction(const Iterate& current_iterate, const Direction& direction) const {
      constexpr double macheps = std::numeric_limits<double>::epsilon();
      for (size_t variable_index: Range(current_iterate.number_variables)) {
         if (std::abs(direction.primals[variable_index]) / (1. + std::abs(current_iterate.primals[variable_index])) >= 10.*macheps) {
            return false;
         }
      }
      if (current_iterate.primal_infeasibility > this->theta_min) {
         return false;
      }
      return true;
   }
} // namespace