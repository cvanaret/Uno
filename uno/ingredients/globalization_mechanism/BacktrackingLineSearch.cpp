// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <cassert>
#include "linear_algebra/SymmetricIndefiniteLinearSystem.hpp"
#include "BacktrackingLineSearch.hpp"
#include "tools/Logger.hpp"

BacktrackingLineSearch::BacktrackingLineSearch(ConstraintRelaxationStrategy& constraint_relaxation_strategy, const Options& options):
      GlobalizationMechanism(constraint_relaxation_strategy),
      backtracking_ratio(options.get_double("LS_backtracking_ratio")),
      minimum_step_length(options.get_double("LS_min_step_length")),
      tolerance(options.get_double("tolerance")),
      use_second_order_correction(options.get_bool("use_second_order_correction")),
      statistics_minor_column_order(options.get_int("statistics_minor_column_order")),
      statistics_SOC_column_order(options.get_int("statistics_SOC_column_order")),
      statistics_LS_step_length_column_order(options.get_int("statistics_LS_step_length_column_order")) {
   // check the initial and minimal step lengths
   assert(0 < this->backtracking_ratio && this->backtracking_ratio < 1. && "The LS backtracking ratio should be in (0, 1)");
   assert(0 < this->minimum_step_length && this->minimum_step_length < 1. && "The LS minimum step length should be in (0, 1)");
}

void BacktrackingLineSearch::initialize(Statistics& statistics, Iterate& first_iterate) {
   statistics.add_column("LS iters", Statistics::int_width + 3, this->statistics_minor_column_order);
   if (this->use_second_order_correction) {
      statistics.add_column("SOC", Statistics::char_width - 2, this->statistics_SOC_column_order);
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
      DEBUG << "Unstable regularization: switching to solving the feasibility problem\n";
      return this->constraint_relaxation_strategy.solve_feasibility_problem(statistics, current_iterate);
   }
}

std::tuple<Iterate, double> BacktrackingLineSearch::compute_acceptable_iterate(Statistics& statistics, const Model& /*model*/, Iterate& current_iterate) {
   // compute the direction
   Direction direction = this->compute_direction(statistics, current_iterate);
   this->solving_feasibility_problem = false;
   this->total_number_iterations = 0;

   // backtrack along the direction
   try {
      return this->backtrack_along_direction(statistics, current_iterate, direction);
   }
   catch (const StepLengthTooSmall& e) {
      // if step length is too small, revert to solving the feasibility problem (if we aren't already solving it)
      if (not this->solving_feasibility_problem && 0. < direction.multipliers.objective && this->tolerance < current_iterate.progress.infeasibility) {
         DEBUG << "The line search terminated with a step length smaller than " << this->minimum_step_length << '\n';
         // reset the line search with the restoration solution
         direction = this->constraint_relaxation_strategy.solve_feasibility_problem(statistics, current_iterate, direction.primals);
         this->solving_feasibility_problem = true;
         auto [trial_iterate, step_norm] = this->backtrack_along_direction(statistics, current_iterate, direction);
         this->total_number_iterations += this->number_iterations;
         return {trial_iterate, step_norm};
      }
      else {
         WARNING << "The feasibility problem failed to make progress\n";
         throw std::runtime_error("Line search: maximum number of iterations reached\n");
      }
   }
}

std::tuple<Iterate, double> BacktrackingLineSearch::backtrack_along_direction(Statistics& statistics, Iterate& current_iterate,
      const Direction& direction) {
   // backtrack on the primal-dual step length computed by the subproblem
   // most subproblem methods return a step length of 1. Interior-point methods however apply the fraction-to-boundary condition
   double primal_dual_step_length = direction.primal_dual_step_length;
   this->number_iterations = 0;

   while (not this->termination(primal_dual_step_length)) {
      this->number_iterations++;
      this->print_iteration(primal_dual_step_length);

      // assemble the trial iterate by going a fraction along the direction
      Iterate trial_iterate = GlobalizationMechanism::assemble_trial_iterate(current_iterate, direction, primal_dual_step_length,
            direction.bound_dual_step_length);
      try {
         if (this->constraint_relaxation_strategy.is_iterate_acceptable(statistics, current_iterate, trial_iterate, direction,
               primal_dual_step_length)) {
            this->total_number_iterations += this->number_iterations;
            this->set_statistics(statistics, direction, primal_dual_step_length);

            double step_norm = primal_dual_step_length * direction.norm;
            return std::make_tuple(std::move(trial_iterate), step_norm);
         }
         // (optional) second-order correction
         else if (this->use_second_order_correction && this->number_iterations == 1 && not this->solving_feasibility_problem &&
                  trial_iterate.progress.infeasibility >= current_iterate.progress.infeasibility) {
            // if full step is rejected, compute a (temporary) SOC direction
            Direction direction_soc = this->constraint_relaxation_strategy.compute_second_order_correction(trial_iterate, direction.primal_dual_step_length);
            if (direction_soc.status != SubproblemStatus::INFEASIBLE) {
               // assemble the (temporary) SOC trial iterate
               Iterate trial_iterate_soc = GlobalizationMechanism::assemble_trial_iterate(current_iterate, direction_soc,
                     direction_soc.primal_dual_step_length, direction.bound_dual_step_length);
               if (this->constraint_relaxation_strategy.is_iterate_acceptable(statistics, current_iterate, trial_iterate_soc, direction_soc,
                     direction_soc.primal_dual_step_length)) {
                  DEBUG << "Trial SOC step accepted\n";
                  this->total_number_iterations += this->number_iterations;
                  this->set_statistics(statistics, direction_soc, direction_soc.primal_dual_step_length);
                  statistics.add_statistic("SOC", "x");

                  // let the subproblem know the accepted iterate
                  trial_iterate_soc.multipliers.lower_bounds = trial_iterate.multipliers.lower_bounds;
                  trial_iterate_soc.multipliers.upper_bounds = trial_iterate.multipliers.upper_bounds;
                  return std::make_tuple(std::move(trial_iterate_soc), direction_soc.norm);
               }
            }
            DEBUG << "Trial SOC step discarded\n\n";
            statistics.add_statistic("SOC", "-");
            this->decrease_step_length(primal_dual_step_length);
         }
         else { // trial iterate not acceptable
            this->decrease_step_length(primal_dual_step_length);
         }
      }
      catch (const EvaluationError& e) {
         GlobalizationMechanism::print_warning(e.what());
         this->decrease_step_length(primal_dual_step_length);
      }
   }
   throw StepLengthTooSmall();
}

void BacktrackingLineSearch::decrease_step_length(double& primal_dual_step_length) const {
   // step length follows the following sequence: 1, ratio, ratio^2, ratio^3, ...
   primal_dual_step_length *= this->backtracking_ratio;
   assert(0 < primal_dual_step_length && primal_dual_step_length <= 1 && "The line-search step length is not in (0, 1]");
}

bool BacktrackingLineSearch::termination(double primal_dual_step_length) const {
   return (primal_dual_step_length < this->minimum_step_length);
}

void BacktrackingLineSearch::set_statistics(Statistics& statistics, const Direction& direction, double primal_dual_step_length) const {
   statistics.add_statistic("LS iters", this->total_number_iterations);
   statistics.add_statistic("LS step length", primal_dual_step_length);
   statistics.add_statistic("step norm", primal_dual_step_length * direction.norm);
}

void BacktrackingLineSearch::print_iteration(double primal_dual_step_length) {
   DEBUG << "\tLINE SEARCH iteration " << this->number_iterations << ", step_length " << primal_dual_step_length << '\n';
}
