// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "BacktrackingLineSearch.hpp"
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategy.hpp"
#include "model/Model.hpp"
#include "optimization/Direction.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "optimization/Iterate.hpp"
#include "ingredients/subproblem_solvers/SubproblemStatus.hpp"
#include "tools/Logger.hpp"
#include "options/Options.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   BacktrackingLineSearch::BacktrackingLineSearch(const Model& model, const Options& options):
         GlobalizationMechanism(model, false, options),
         backtracking_ratio(options.get_double("LS_backtracking_ratio")),
         minimum_step_length(options.get_double("LS_min_step_length")),
         scale_duals_with_step_length(options.get_bool("LS_scale_duals_with_step_length")) {
      // check the initial and minimal step lengths
      assert(0 < this->backtracking_ratio && this->backtracking_ratio < 1. && "The LS backtracking ratio should be in (0, 1)");
      assert(0 < this->minimum_step_length && this->minimum_step_length < 1. && "The LS minimum step length should be in (0, 1)");
   }

   void BacktrackingLineSearch::initialize(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Direction& direction) {
      this->constraint_relaxation_strategy->initialize(statistics, model, current_iterate, direction, INF<double>);
      statistics.add_column("Minor", Statistics::int_width, 3, Statistics::column_order.at("Minor"));
      statistics.add_column("Steplength", Statistics::double_width + 1, 2, Statistics::column_order.at("Steplength"));
   }

   void BacktrackingLineSearch::compute_next_iterate(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Iterate& trial_iterate, Direction& direction, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) {
      DEBUG2 << "Current iterate\n" << current_iterate << '\n';

      this->constraint_relaxation_strategy->compute_feasible_direction(statistics, current_iterate, direction, INF<double>,
         warmstart_information);
      BacktrackingLineSearch::check_unboundedness(direction);
      this->backtrack_along_direction(statistics, model, current_iterate, trial_iterate, direction, warmstart_information,
         user_callbacks);
   }

   std::string BacktrackingLineSearch::get_name() const {
      return "LS " + this->constraint_relaxation_strategy->get_name();
   }

   // protected member functions

   // go a fraction along the direction by finding an acceptable step length
   void BacktrackingLineSearch::backtrack_along_direction(Statistics& statistics, const Model& model, Iterate& current_iterate,
         Iterate& trial_iterate, Direction& direction, WarmstartInformation& warmstart_information, UserCallbacks& user_callbacks) const {
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
            GlobalizationMechanism::assemble_trial_iterate(model, current_iterate, trial_iterate, direction, step_length,
               // scale or not the constraint dual direction with the LS step length
               this->scale_duals_with_step_length ? step_length : 1.);
            statistics.set("||Step||", step_length * direction.norm);

            is_acceptable = this->constraint_relaxation_strategy->is_iterate_acceptable(statistics, model, current_iterate,
               trial_iterate, direction, step_length, warmstart_information, user_callbacks);
            GlobalizationMechanism::set_primal_statistics(statistics, model, trial_iterate);
         }
         catch (const EvaluationError&) {
            statistics.set("Status", "eval. error");
         }
         BacktrackingLineSearch::set_LS_statistics(statistics, number_iterations);

         if (is_acceptable) {
            trial_iterate.status = this->constraint_relaxation_strategy->check_termination(model, trial_iterate);
            GlobalizationMechanism::set_dual_residuals_statistics(statistics, trial_iterate);
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
            termination = BacktrackingLineSearch::terminate_with_small_step_length(statistics, model, trial_iterate);
            if (!termination) {
               // test if we can switch to solving the feasibility problem
               if (this->constraint_relaxation_strategy->solving_feasibility_problem() || !model.is_constrained()) {
                  throw std::runtime_error("LS failed");
               }
               // switch to solving the feasibility problem
               statistics.set("Status", "small step length");
               this->constraint_relaxation_strategy->switch_to_feasibility_problem(statistics, current_iterate, INF<double>,
                  warmstart_information);
               this->constraint_relaxation_strategy->compute_feasible_direction(statistics, current_iterate, direction,
                  INF<double>, warmstart_information);
               BacktrackingLineSearch::check_unboundedness(direction);
               // restart backtracking
               step_length = 1.;
               number_iterations = 0;
            }
         }
      } // end while loop
   }

   bool BacktrackingLineSearch::terminate_with_small_step_length(Statistics& statistics, const Model& model,
         Iterate& trial_iterate) const {
      bool termination = false;
      trial_iterate.status = this->constraint_relaxation_strategy->check_termination(model, trial_iterate);
      if (trial_iterate.status != SolutionStatus::NOT_OPTIMAL) {
         statistics.set("Status", "accepted (small step length)");
         GlobalizationMechanism::set_dual_residuals_statistics(statistics, trial_iterate);
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
         throw std::runtime_error("The subproblem is unbounded, this should not happen. If the subproblem has curvature,"
            "use regularization. If not, use a trust-region method.\n");
      }
   }

   void BacktrackingLineSearch::set_LS_statistics(Statistics& statistics, size_t number_iterations) {
      statistics.set("Minor", number_iterations);
   }
} // namespace