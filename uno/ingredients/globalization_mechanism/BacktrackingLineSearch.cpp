// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <cassert>
#include "linear_algebra/SymmetricIndefiniteLinearSystem.hpp"
#include "BacktrackingLineSearch.hpp"
#include "tools/Logger.hpp"

BacktrackingLineSearch::BacktrackingLineSearch(ConstraintRelaxationStrategy& constraint_relaxation_strategy, const Options& options):
      GlobalizationMechanism(constraint_relaxation_strategy),
      backtracking_ratio(options.get_double("LS_backtracking_ratio")),
      min_step_length(options.get_double("LS_min_step_length")),
      use_second_order_correction(options.get_bool("use_second_order_correction")),
      statistics_SOC_column_order(options.get_int("statistics_SOC_column_order")),
      statistics_LS_step_length_column_order(options.get_int("statistics_LS_step_length_column_order")) {
   // check the initial and minimal step lengths
   assert(0 < this->backtracking_ratio && this->backtracking_ratio < 1. && "The LS backtracking ratio should be in (0, 1)");
   assert(0 < this->min_step_length && this->min_step_length < 1. && "The LS minimum step length should be in (0, 1)");
}

void BacktrackingLineSearch::initialize(Statistics& statistics, Iterate& first_iterate) {
   if (this->use_second_order_correction) {
      statistics.add_column("SOC", Statistics::char_width, this->statistics_SOC_column_order);
   }
   statistics.add_column("LS step length", Statistics::double_width, this->statistics_LS_step_length_column_order);

   // generate the initial point
   this->constraint_relaxation_strategy.initialize(statistics, first_iterate);
}

Direction BacktrackingLineSearch::compute_direction(Statistics& statistics, Iterate& current_iterate) {
   try {
      this->solving_feasibility_problem = false;
      return this->constraint_relaxation_strategy.compute_feasible_direction(statistics, current_iterate);
   }
   catch (const UnstableRegularization&) {
      this->solving_feasibility_problem = true;
      return this->constraint_relaxation_strategy.solve_feasibility_problem(statistics, current_iterate);
   }
}

std::tuple<Iterate, double> BacktrackingLineSearch::compute_acceptable_iterate(Statistics& statistics, const Model& /*model*/, Iterate& current_iterate) {
   // compute the direction
   Direction direction = this->compute_direction(statistics, current_iterate);
   this->solving_feasibility_problem = false;

   // backtrack along the direction
   this->step_length = 1.;
   this->number_iterations = 0;
   this->total_number_iterations = 0;
   bool failure = false;
   while (not failure) {
      while (not this->termination()) {
         this->number_iterations++;
         this->print_iteration();

         // assemble the trial iterate by going a fraction along the direction
         Iterate trial_iterate = GlobalizationMechanism::assemble_trial_iterate(current_iterate, direction, this->step_length);
         try {
            const bool is_acceptable = this->constraint_relaxation_strategy.is_iterate_acceptable(statistics, current_iterate, trial_iterate,
                  direction, this->step_length);
            if (is_acceptable) {
               this->total_number_iterations += this->number_iterations;
               this->constraint_relaxation_strategy.postprocess_accepted_iterate(trial_iterate);
               this->set_statistics(statistics, direction);

               double step_norm = this->step_length * direction.norm;
               return std::make_tuple(std::move(trial_iterate), step_norm);
            }
            // (optional) second-order correction
            else if (this->use_second_order_correction && this->number_iterations == 1 && not this->solving_feasibility_problem &&
                  trial_iterate.nonlinear_progress.infeasibility >= current_iterate.nonlinear_progress.infeasibility) {
               // if full step is rejected, compute a (temporary) SOC direction
               Direction direction_soc = this->constraint_relaxation_strategy.compute_second_order_correction(trial_iterate);
               if (direction_soc.status == SubproblemStatus::INFEASIBLE) {
                  DEBUG << "Trial SOC subproblem infeasible\n\n";
                  statistics.add_statistic("SOC", "-");
                  this->decrease_step_length();
               }
               else {
                  // assemble the (temporary) SOC trial iterate
                  Iterate trial_iterate_soc = GlobalizationMechanism::assemble_trial_iterate(current_iterate, direction_soc, this->step_length);
                  const bool is_SOC_acceptable = this->constraint_relaxation_strategy.is_iterate_acceptable(statistics, current_iterate,
                        trial_iterate_soc, direction_soc, this->step_length);
                  if (is_SOC_acceptable) {
                     DEBUG << "Trial SOC step accepted\n";
                     this->set_statistics(statistics, direction_soc);
                     statistics.add_statistic("SOC", "x");

                     // let the subproblem know the accepted iterate
                     this->constraint_relaxation_strategy.postprocess_accepted_iterate(trial_iterate_soc);
                     trial_iterate_soc.multipliers.lower_bounds = trial_iterate.multipliers.lower_bounds;
                     trial_iterate_soc.multipliers.upper_bounds = trial_iterate.multipliers.upper_bounds;
                     return std::make_tuple(std::move(trial_iterate_soc), direction_soc.norm);
                  }
                  else {
                     DEBUG << "Trial SOC step discarded\n\n";
                     statistics.add_statistic("SOC", "-");
                     this->decrease_step_length();
                  }
               }
            }
            else { // trial iterate not acceptable
               this->decrease_step_length();
            }
         }
         catch (const EvaluationError& e) {
            GlobalizationMechanism::print_warning(e.what());
            this->decrease_step_length();
         }
      }
      // if step length is too small, revert to solving the feasibility problem (if we aren't already solving it)
      if (not this->solving_feasibility_problem && 0. < direction.multipliers.objective) {
         // TODO: test if 0. < current_iterate.progress.infeasibility ?
         DEBUG << "The line search failed, switching to feasibility problem\n";
         // reset the line search with the restoration solution
         direction = this->constraint_relaxation_strategy.solve_feasibility_problem(statistics, current_iterate, direction.primals);
         this->step_length = 1.;
         this->total_number_iterations += this->number_iterations;
         this->number_iterations = 0;
         this->solving_feasibility_problem = true;
      }
      else {
         WARNING << "The feasibility problem failed to make progress\n";
         failure = true;
      }
   }
   throw std::runtime_error("Line search: maximum number of iterations reached\n");
}

void BacktrackingLineSearch::decrease_step_length() {
   // step length follows the following sequence: 1, ratio, ratio^2, ratio^3, ...
   this->step_length *= this->backtracking_ratio;
   assert(0 < this->step_length && this->step_length <= 1 && "The line-search step length is not in (0, 1]");
}

bool BacktrackingLineSearch::termination() const {
   return (this->step_length < this->min_step_length);
}

void BacktrackingLineSearch::set_statistics(Statistics& statistics, const Direction& direction) const {
   statistics.add_statistic("minor", this->total_number_iterations);
   statistics.add_statistic("LS step length", this->step_length);
   statistics.add_statistic("step norm", this->step_length * direction.norm);
}

void BacktrackingLineSearch::print_iteration() {
   DEBUG << "\tLINE SEARCH iteration " << this->number_iterations << ", step_length " << this->step_length << '\n';
}
