// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategy.hpp"
#include "BacktrackingLineSearch.hpp"
#include "model/Model.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "tools/Logger.hpp"
#include "options/Options.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   BacktrackingLineSearch::BacktrackingLineSearch(ConstraintRelaxationStrategy& constraint_relaxation_strategy, const Options& options):
         GlobalizationMechanism(constraint_relaxation_strategy),
         backtracking_ratio(options.get_double("LS_backtracking_ratio")),
         minimum_step_length(options.get_double("LS_min_step_length")),
         scale_duals_with_step_length(options.get_bool("LS_scale_duals_with_step_length")) {
      // check the initial and minimal step lengths
      assert(0 < this->backtracking_ratio && this->backtracking_ratio < 1. && "The LS backtracking ratio should be in (0, 1)");
      assert(0 < this->minimum_step_length && this->minimum_step_length < 1. && "The LS minimum step length should be in (0, 1)");
   }

   void BacktrackingLineSearch::initialize(Statistics& statistics, Iterate& initial_iterate, const Options& options) {
      statistics.add_column("LS iter", Statistics::int_width + 2, options.get_int("statistics_minor_column_order"));
      statistics.add_column("step length", Statistics::double_width - 4, options.get_int("statistics_LS_step_length_column_order"));
      
      this->constraint_relaxation_strategy.initialize(statistics, initial_iterate, options);
   }

   void BacktrackingLineSearch::compute_next_iterate(Statistics& statistics, const Model& model, Iterate& current_iterate, Iterate& trial_iterate,
         WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) {
      DEBUG2 << "Current iterate\n" << current_iterate << '\n';

      this->constraint_relaxation_strategy.compute_feasible_direction(statistics, current_iterate, this->direction, warmstart_information);
      BacktrackingLineSearch::check_unboundedness(this->direction);
      this->backtrack_along_direction(statistics, model, current_iterate, trial_iterate, warmstart_information, user_callbacks);
   }

   // go a fraction along the direction by finding an acceptable step length
   void BacktrackingLineSearch::backtrack_along_direction(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Iterate& trial_iterate, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) {
      double step_length = 1.;
      bool termination = false;
      size_t number_iterations = 0;
      while (not termination) {
         number_iterations++;
         DEBUG << "\n\tLine-search iteration " << number_iterations << ", step_length " << step_length << '\n';
         if (1 < number_iterations) { statistics.start_new_line(); }
         statistics.set("step length", step_length);

         bool is_acceptable = false;
         try {
            // take a step as a fraction of the direction
            GlobalizationMechanism::assemble_trial_iterate(model, current_iterate, trial_iterate, this->direction, step_length,
                  // scale or not the constraint dual direction with the LS step length
                  this->scale_duals_with_step_length ? step_length : 1.);

            is_acceptable = this->constraint_relaxation_strategy.is_iterate_acceptable(statistics, current_iterate, trial_iterate, this->direction,
                  step_length, warmstart_information, user_callbacks);
            this->set_statistics(statistics, trial_iterate, this->direction, step_length, number_iterations);
         }
         catch (const EvaluationError& e) {
            this->set_statistics(statistics, number_iterations);
            statistics.set("status", "eval. error");
         }

         if (is_acceptable) {
            trial_iterate.status = this->constraint_relaxation_strategy.check_termination(trial_iterate);
            this->constraint_relaxation_strategy.set_dual_residuals_statistics(statistics, trial_iterate);
            termination = true;
            if (Logger::level == INFO) statistics.print_current_line();
         }
         else if (step_length >= this->minimum_step_length) {
            step_length = this->decrease_step_length(step_length);
            if (Logger::level == INFO) statistics.print_current_line();
         }
         else { // minimum_step_length reached
            DEBUG << "The line search step length is smaller than " << this->minimum_step_length << '\n';
            // check if we can terminate at a first-order point
            termination = this->terminate_with_small_step_length(statistics, trial_iterate);
            if (not termination) {
               // test if we can switch to solving the feasibility problem
               if (this->constraint_relaxation_strategy.solving_feasibility_problem() || not model.is_constrained()) {
                  throw std::runtime_error("LS failed");
               }
               // switch to solving the feasibility problem
               statistics.set("status", "small step length");
               this->constraint_relaxation_strategy.switch_to_feasibility_problem(statistics, current_iterate, warmstart_information);
               this->constraint_relaxation_strategy.compute_feasible_direction(statistics, current_iterate, this->direction, this->direction.primals,
                     warmstart_information);
               BacktrackingLineSearch::check_unboundedness(this->direction);
               // restart backtracking
               step_length = 1.;
               number_iterations = 0;
            }
         }
      } // end while loop
   }

   bool BacktrackingLineSearch::terminate_with_small_step_length(Statistics& statistics, Iterate& trial_iterate) {
      bool termination = false;
      trial_iterate.status = this->constraint_relaxation_strategy.check_termination(trial_iterate);
      if (trial_iterate.status != IterateStatus::NOT_OPTIMAL) {
         statistics.set("status", "accepted (small step length)");
         this->constraint_relaxation_strategy.set_dual_residuals_statistics(statistics, trial_iterate);
         termination = true;
      }
      return termination;
   }

   // step length follows the following sequence: 1, ratio, ratio^2, ratio^3, ...
   double BacktrackingLineSearch::decrease_step_length(double step_length) const {
      step_length *= this->backtracking_ratio;
      assert(0 < step_length && step_length <= 1 && "The line-search step length is not in (0, 1]");
      return step_length;
   }

   void BacktrackingLineSearch::check_unboundedness(const Direction& direction) {
      if (direction.status == SubproblemStatus::UNBOUNDED_PROBLEM) {
         throw std::runtime_error("The subproblem is unbounded, this should not happen. If the subproblem has curvature, use regularization. If not, "
                                  "use a trust-region method.\n");
      }
   }

   void BacktrackingLineSearch::set_statistics(Statistics& statistics, size_t number_iterations) const {
      statistics.set("LS iter", number_iterations);
   }

   void BacktrackingLineSearch::set_statistics(Statistics& statistics, const Iterate& trial_iterate, const Direction& direction,
         double primal_dual_step_length, size_t number_iterations) const {
      if (trial_iterate.is_objective_computed) {
         statistics.set("objective", trial_iterate.evaluations.objective);
      }
      statistics.set("step norm", primal_dual_step_length * direction.norm);
      this->set_statistics(statistics, number_iterations);
   }
} // namespace
