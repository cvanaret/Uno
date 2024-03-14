// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "linear_algebra/SymmetricIndefiniteLinearSystem.hpp"
#include "BacktrackingLineSearch.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "tools/Logger.hpp"

BacktrackingLineSearch::BacktrackingLineSearch(Statistics& statistics, ConstraintRelaxationStrategy& constraint_relaxation_strategy,
         const Options& options):
      GlobalizationMechanism(constraint_relaxation_strategy, options),
      backtracking_ratio(options.get_double("LS_backtracking_ratio")),
      minimum_step_length(options.get_double("LS_min_step_length")),
      scale_duals_with_step_length(options.get_bool("LS_scale_duals_with_step_length")) {
   // check the initial and minimal step lengths
   assert(0 < this->backtracking_ratio && this->backtracking_ratio < 1. && "The LS backtracking ratio should be in (0, 1)");
   assert(0 < this->minimum_step_length && this->minimum_step_length < 1. && "The LS minimum step length should be in (0, 1)");

   statistics.add_column("LS iter", Statistics::int_width + 2, options.get_int("statistics_minor_column_order"));
   statistics.add_column("step length", Statistics::double_width - 3, options.get_int("statistics_LS_step_length_column_order"));
}

void BacktrackingLineSearch::initialize(Statistics& statistics, Iterate& initial_iterate) {
   this->constraint_relaxation_strategy.initialize(statistics, initial_iterate);
}

Iterate BacktrackingLineSearch::compute_next_iterate(Statistics& statistics, const Model& model, Iterate& current_iterate) {
   WarmstartInformation warmstart_information{};
   warmstart_information.set_hot_start();
   DEBUG2 << "Current iterate\n" << current_iterate << '\n';

   // compute the direction
   Direction direction = this->constraint_relaxation_strategy.compute_feasible_direction(statistics, current_iterate, warmstart_information);
   BacktrackingLineSearch::check_unboundedness(direction);

   // backtrack along the direction
   this->total_number_iterations = 0;
   return this->backtrack_along_direction(statistics, model, current_iterate, direction, warmstart_information);
}

// backtrack on the primal-dual step length computed by the subproblem
Iterate BacktrackingLineSearch::backtrack_along_direction(Statistics& statistics, const Model& model, Iterate& current_iterate,
      const Direction& direction, WarmstartInformation& warmstart_information) {
   // most subproblem methods return a step length of 1. Interior-point methods however apply the fraction-to-boundary condition
   double step_length = direction.primal_dual_step_length;
   bool reached_small_step_length = false;
   size_t number_iterations = 0;
   while (not reached_small_step_length) {
      this->total_number_iterations++;
      number_iterations++;
      if (1 < number_iterations) {
         statistics.new_line();
      }
      BacktrackingLineSearch::print_iteration(number_iterations, step_length);
      statistics.add_statistic("step length", step_length);

      try {
         // assemble the trial iterate by going a fraction along the direction
         Iterate trial_iterate = this->assemble_trial_iterate(model, current_iterate, direction, step_length);
         // check whether the trial iterate is accepted
         if (this->constraint_relaxation_strategy.is_iterate_acceptable(statistics, current_iterate, trial_iterate, direction, step_length)) {
            // check termination criteria
            trial_iterate.status = this->check_convergence(model, trial_iterate);
            this->set_statistics(statistics, direction, step_length);
            if (Logger::level == INFO) statistics.print_current_line();
            return trial_iterate;
         }
         // small step length
         else if (step_length < this->minimum_step_length) {
            DEBUG << "The line search step length is smaller than " << this->minimum_step_length << '\n';
            reached_small_step_length = true;
            this->set_statistics(statistics, direction, step_length);
         }
         else {
            this->set_statistics(statistics, direction, step_length);
            step_length = this->decrease_step_length(step_length);
         }
         if (Logger::level == INFO) statistics.print_current_line();
      }
      catch (const EvaluationError& e) {
         WARNING << YELLOW << e.what() << RESET;
         this->set_statistics(statistics);
         statistics.add_statistic("status", "eval. error");
         if (Logger::level == INFO) statistics.print_current_line();
         step_length = this->decrease_step_length(step_length);
      }
   }

   // reached a small step length: revert to solving the feasibility problem
   if (this->constraint_relaxation_strategy.solving_feasibility_problem()) {
      throw std::runtime_error("LS on feasibility problem failed.");
   }
   else {
      warmstart_information.set_cold_start();
      this->constraint_relaxation_strategy.switch_to_feasibility_problem(statistics, current_iterate, warmstart_information);
      Direction direction_feasibility = this->constraint_relaxation_strategy.compute_feasible_direction(statistics, current_iterate,
            direction.primals, warmstart_information);
      BacktrackingLineSearch::check_unboundedness(direction_feasibility);
      Iterate trial_iterate_feasibility = this->backtrack_along_direction(statistics, model, current_iterate, direction_feasibility,
            warmstart_information);
      return trial_iterate_feasibility;
   }
}

Iterate BacktrackingLineSearch::assemble_trial_iterate(const Model& model, Iterate& current_iterate, const Direction& direction,
      double primal_dual_step_length) const {
   Iterate trial_iterate = GlobalizationMechanism::assemble_trial_iterate(current_iterate, direction, primal_dual_step_length,
         // scale or not the dual directions with the step lengths
         this->scale_duals_with_step_length ? primal_dual_step_length : 1.,
         this->scale_duals_with_step_length ? direction.bound_dual_step_length : 1.);

   // project the steps within the bounds to avoid numerical errors
   model.project_primals_onto_bounds(trial_iterate.primals);
   return trial_iterate;
}

double BacktrackingLineSearch::decrease_step_length(double step_length) const {
   // step length follows the following sequence: 1, ratio, ratio^2, ratio^3, ...
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
   statistics.add_statistic("LS iter", this->total_number_iterations);
}

void BacktrackingLineSearch::set_statistics(Statistics& statistics, const Direction& direction, double primal_dual_step_length) const {
   statistics.add_statistic("step norm", primal_dual_step_length * direction.norm);
   this->set_statistics(statistics);
}

void BacktrackingLineSearch::print_iteration(size_t number_iterations, double primal_dual_step_length) {
   DEBUG << "\tLINE SEARCH iteration " << number_iterations << ", step_length " << primal_dual_step_length << '\n';
}
