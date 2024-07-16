// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "ingredients/constraint_relaxation_strategy/ConstraintRelaxationStrategy.hpp"
#include "BacktrackingLineSearch.hpp"
#include "model/Model.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "tools/Logger.hpp"
#include "tools/Options.hpp"
#include "tools/Statistics.hpp"

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

void BacktrackingLineSearch::compute_next_iterate(Statistics& statistics, const Model& model, Iterate& current_iterate, Iterate& trial_iterate) {
   WarmstartInformation warmstart_information{};
   warmstart_information.set_hot_start();
   DEBUG2 << "Current iterate\n" << current_iterate << '\n';

   // compute the direction
   this->constraint_relaxation_strategy.compute_feasible_direction(statistics, current_iterate, this->direction, warmstart_information);
   BacktrackingLineSearch::check_unboundedness(this->direction);

   // backtrack along the direction
   this->total_number_iterations = 0;
   this->backtrack_along_direction(statistics, model, current_iterate, trial_iterate, this->direction, warmstart_information);
}

// backtrack on the primal-dual step length computed by the subproblem
void BacktrackingLineSearch::backtrack_along_direction(Statistics& statistics, const Model& model, Iterate& current_iterate,
      Iterate& trial_iterate, const Direction& direction, WarmstartInformation& warmstart_information) {
   // most subproblem methods return a step length of 1. Interior-point methods however apply the fraction-to-boundary condition
   double step_length = 1.;
   bool reached_small_step_length = false;
   size_t number_iterations = 0;
   while (not reached_small_step_length) {
      this->total_number_iterations++;
      number_iterations++;
      if (1 < number_iterations) {
         statistics.start_new_line();
      }
      BacktrackingLineSearch::print_iteration(number_iterations, step_length);
      statistics.set("step length", step_length);

      try {
         // assemble the trial iterate by going a fraction along the direction
         GlobalizationMechanism::assemble_trial_iterate(model, current_iterate, trial_iterate, direction, step_length,
               // scale or not the constraint dual direction with the LS step length
               this->scale_duals_with_step_length ? step_length : 1.);

         // check whether the trial iterate is accepted
         if (this->constraint_relaxation_strategy.is_iterate_acceptable(statistics, current_iterate, trial_iterate, direction, step_length)) {
            this->set_statistics(statistics, trial_iterate, direction, step_length);
            if (Logger::level == INFO) statistics.print_current_line();
            return;
         }
         // small step length
         else if (step_length < this->minimum_step_length) {
            DEBUG << "The line search step length is smaller than " << this->minimum_step_length << '\n';
            reached_small_step_length = true;
            this->set_statistics(statistics, trial_iterate, direction, step_length);
         }
         else {
            this->set_statistics(statistics, trial_iterate, direction, step_length);
            step_length = this->decrease_step_length(step_length);
         }
         if (Logger::level == INFO) statistics.print_current_line();
      }
      catch (const EvaluationError& e) {
         this->set_statistics(statistics);
         statistics.set("status", "eval. error");
         if (Logger::level == INFO) statistics.print_current_line();
         step_length = this->decrease_step_length(step_length);
      }
   }

   // reached a small step length: revert to solving the feasibility problem
   if (this->constraint_relaxation_strategy.solving_feasibility_problem()) {
      throw std::runtime_error("Feasibility LS failed");
   }
   else if (not model.is_constrained()) {
      throw std::runtime_error("Regular LS failed");
   }
   else {
      warmstart_information.set_cold_start();
      statistics.set("status", "LS failed");
      this->constraint_relaxation_strategy.switch_to_feasibility_problem(statistics, current_iterate);
      warmstart_information.set_cold_start();
      this->constraint_relaxation_strategy.compute_feasible_direction(statistics, current_iterate, this->direction, this->direction.primals,
            warmstart_information);
      BacktrackingLineSearch::check_unboundedness(this->direction);
      this->backtrack_along_direction(statistics, model, current_iterate, trial_iterate, this->direction, warmstart_information);
   }
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

void BacktrackingLineSearch::set_statistics(Statistics& statistics) const {
   statistics.set("LS iter", this->total_number_iterations);
}

void BacktrackingLineSearch::set_statistics(Statistics& statistics, const Iterate& trial_iterate, const Direction& direction,
      double primal_dual_step_length) const {
   if (trial_iterate.is_objective_computed) {
      statistics.set("objective", trial_iterate.evaluations.objective);
   }
   statistics.set("step norm", primal_dual_step_length * direction.norm);
   this->set_statistics(statistics);
}

void BacktrackingLineSearch::print_iteration(size_t number_iterations, double primal_dual_step_length) {
   DEBUG << "\n\tLINE SEARCH iteration " << number_iterations << ", step_length " << primal_dual_step_length << '\n';
}
