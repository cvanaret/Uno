// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "BacktrackingLineSearch.hpp"
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategy.hpp"
#include "ingredients/inertia_correction_strategies/UnstableInertiaCorrection.hpp"
#include "model/Model.hpp"
#include "optimization/Direction.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "optimization/Iterate.hpp"
#include "ingredients/subproblem_solvers/SubproblemStatus.hpp"
#include "optimization/EvaluationCache.hpp"
#include "tools/Logger.hpp"
#include "options/Options.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   BacktrackingLineSearch::BacktrackingLineSearch(const Model& model, Options& options):
         GlobalizationMechanism(model, false, options),
         backtracking_ratio(options.get_double("LS_backtracking_ratio")),
         minimum_step_length(options.get_double("LS_min_step_length")),
         scale_duals_with_step_length(options.get_bool("LS_scale_duals_with_step_length")) {
      // check the initial and minimal step lengths
      assert(0 < this->backtracking_ratio && this->backtracking_ratio < 1. && "The LS backtracking ratio should be in (0, 1)");
      assert(0 < this->minimum_step_length && this->minimum_step_length < 1. && "The LS minimum step length should be in (0, 1)");
   }

   void BacktrackingLineSearch::initialize(Statistics& statistics, const Model& model, Iterate& current_iterate,
         EvaluationCache& evaluation_cache, Options& options) {
      this->constraint_relaxation_strategy->initialize(statistics, model, current_iterate, this->direction, false,
         evaluation_cache, options);
      statistics.add_column("Minor", Statistics::int_width, 3, Statistics::column_order.at("Minor"));
      statistics.add_column("Steplength", Statistics::double_width + 1, 2, Statistics::column_order.at("Steplength"));
      GlobalizationMechanism::set_primal_statistics(statistics, model, current_iterate, evaluation_cache.current_evaluations);
      GlobalizationMechanism::set_dual_residuals_statistics(statistics, current_iterate);
   }

   void BacktrackingLineSearch::compute_next_iterate(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Iterate& trial_iterate, EvaluationCache& evaluation_cache, WarmstartInformation& warmstart_information,
         UserCallbacks& user_callbacks) {
      DEBUG2 << "Current iterate\n" << current_iterate << '\n';

      // compute a feasible direction
      try {
         this->constraint_relaxation_strategy->compute_feasible_direction(statistics, current_iterate, this->direction,
            INF<double>, evaluation_cache.current_evaluations, warmstart_information);
         BacktrackingLineSearch::check_unboundedness(this->direction);
         const bool backtracking_success = this->backtrack_along_direction(statistics, model, current_iterate, trial_iterate,
            this->direction, evaluation_cache, warmstart_information, user_callbacks);
         if (backtracking_success) {
            return;
         }
         else {
            // if the line search failed, switch to solving the feasibility problem (test first if we can)
            if (this->constraint_relaxation_strategy->solving_feasibility_problem() || !model.is_constrained()) {
               throw std::runtime_error("The line search failed");
            }
            this->constraint_relaxation_strategy->switch_to_feasibility_problem(statistics, current_iterate, this->direction,
               evaluation_cache.current_evaluations, warmstart_information);
         }
      }
      // if the inertia correction failed, switch to solving the feasibility problem
      catch (const UnstableInertiaCorrection&) {
         this->constraint_relaxation_strategy->switch_to_feasibility_problem(statistics, current_iterate, this->direction,
            evaluation_cache.current_evaluations, warmstart_information);
      }

      // solve the feasibility problem
      assert(this->constraint_relaxation_strategy->solving_feasibility_problem());
      this->constraint_relaxation_strategy->compute_feasible_direction(statistics, current_iterate, this->direction,
         INF<double>, evaluation_cache.current_evaluations, warmstart_information);
      BacktrackingLineSearch::check_unboundedness(this->direction);
      const bool backtracking_success = this->backtrack_along_direction(statistics, model, current_iterate,
         trial_iterate, this->direction, evaluation_cache, warmstart_information, user_callbacks);
      if (!backtracking_success) {
         throw std::runtime_error("The line search failed");
      }
   }

   std::string BacktrackingLineSearch::get_name() const {
      return "LS " + this->constraint_relaxation_strategy->get_name();
   }

   // protected member functions

   // go a fraction along the direction by finding an acceptable step length
   // returns true upon success, false upon failure
   bool BacktrackingLineSearch::backtrack_along_direction(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Iterate& trial_iterate, const Direction& direction, EvaluationCache& evaluation_cache, WarmstartInformation& warmstart_information,
         UserCallbacks& user_callbacks) const {
      double step_length = 1.;
      bool termination = false;
      size_t number_iterations = 0;
      while (!termination) {
         ++number_iterations;
         DEBUG << "\n\tLine-search iteration " << number_iterations << ", step_length " << step_length << '\n';
         if (1 < number_iterations) { statistics.start_new_line(); }
         statistics.set("Steplength", step_length);

         bool is_acceptable = false;
         try {
            // take a step as a fraction of the direction
            BacktrackingLineSearch::assemble_trial_iterate(model, current_iterate, trial_iterate, direction, step_length,
               // scale or not the constraint dual direction with the LS step length
               this->scale_duals_with_step_length ? step_length : 1.);
            statistics.set("||Step||", step_length * direction.norm);

            is_acceptable = this->constraint_relaxation_strategy->is_iterate_acceptable(statistics, model, current_iterate,
               trial_iterate, direction, step_length, evaluation_cache, warmstart_information, user_callbacks);
            BacktrackingLineSearch::set_primal_statistics(statistics, model, trial_iterate, evaluation_cache.trial_evaluations);
         }
         catch (const EvaluationError&) {
            statistics.set("Status", "eval. error");
         }
         statistics.set("Minor", number_iterations);

         if (is_acceptable) {
            termination = true;
            GlobalizationMechanism::set_dual_residuals_statistics(statistics, trial_iterate);
            if (Logger::level == INFO) statistics.print_current_line();
         }
         else if (step_length >= this->minimum_step_length) {
            step_length = this->decrease_step_length(step_length);
            evaluation_cache.trial_evaluations.reset();
            if (Logger::level == INFO) statistics.print_current_line();
         }
         else { // minimum_step_length reached
            DEBUG << "The line search step length is smaller than " << this->minimum_step_length << '\n';
            // check if we can terminate at a first-order point
            if (trial_iterate.status != SolutionStatus::NOT_OPTIMAL) {
               statistics.set("Status", "accepted (small step length)");
               termination = true;
            }
            else {
               // switch to solving the feasibility problem
               statistics.set("Status", "small step length");
               return false;
            }
         }
      } // end while loop
      return true;
   }

   // step length follows the following sequence: 1, ratio, ratio^2, ratio^3, ...
   double BacktrackingLineSearch::decrease_step_length(double step_length) const {
      step_length *= this->backtracking_ratio;
      assert(0 < step_length && step_length <= 1 && "The line-search step length is not in (0, 1]");
      return step_length;
   }

   void BacktrackingLineSearch::check_unboundedness(const Direction& direction) {
      if (direction.status == SubproblemStatus::UNBOUNDED_PROBLEM) {
         throw std::runtime_error("The subproblem is unbounded, this should not happen. If the subproblem has curvature,"
            "use regularization. If not, use a trust-region method.\n");
      }
   }
} // namespace