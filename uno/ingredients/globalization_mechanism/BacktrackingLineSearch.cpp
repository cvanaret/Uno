// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
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
      tolerance(options.get_double("tolerance")) {
   // check the initial and minimal step lengths
   assert(0 < this->backtracking_ratio && this->backtracking_ratio < 1. && "The LS backtracking ratio should be in (0, 1)");
   assert(0 < this->minimum_step_length && this->minimum_step_length < 1. && "The LS minimum step length should be in (0, 1)");

   statistics.add_column("LS iters", Statistics::int_width + 3, options.get_int("statistics_minor_column_order"));
   statistics.add_column("LS step length", Statistics::double_width, options.get_int("statistics_LS_step_length_column_order"));
}

void BacktrackingLineSearch::initialize(Iterate& initial_iterate) {
   this->constraint_relaxation_strategy.initialize(initial_iterate);
}

Iterate BacktrackingLineSearch::compute_next_iterate(Statistics& statistics, const Model& model, Iterate& current_iterate) {
   WarmstartInformation warmstart_information{};
   warmstart_information.set_hot_start();
   DEBUG2 << "Current iterate\n" << current_iterate << '\n';

   // compute the direction
   Direction direction = this->compute_direction(statistics, current_iterate, warmstart_information);
   BacktrackingLineSearch::check_unboundedness(direction);
   this->solving_feasibility_problem = false;
   this->total_number_iterations = 0;

   // backtrack along the direction
   return this->backtrack_along_direction(statistics, model, current_iterate, direction, warmstart_information);
}

Direction BacktrackingLineSearch::compute_direction(Statistics& statistics, Iterate& current_iterate, WarmstartInformation& warmstart_information) {
   try {
      this->solving_feasibility_problem = false;
      return this->constraint_relaxation_strategy.compute_feasible_direction(statistics, current_iterate, warmstart_information);
   }
   catch (const UnstableRegularization&) {
      this->solving_feasibility_problem = true;
      DEBUG << "Unstable regularization: switching to solving the feasibility problem\n";
      return this->constraint_relaxation_strategy.solve_feasibility_problem(statistics, current_iterate, warmstart_information);
   }
}

// backtrack on the primal-dual step length computed by the subproblem
Iterate BacktrackingLineSearch::backtrack_along_direction(Statistics& statistics, const Model& model, Iterate& current_iterate,
      const Direction& direction, WarmstartInformation& warmstart_information) {
   // most subproblem methods return a step length of 1. Interior-point methods however apply the fraction-to-boundary condition
   double step_length = direction.primal_dual_step_length;
   this->number_iterations = 0;

   while (not this->termination(step_length)) {
      this->number_iterations++;
      this->print_iteration(step_length);

      try {
         // assemble the trial iterate by going a fraction along the direction
         Iterate trial_iterate = this->assemble_trial_iterate(model, current_iterate, direction, step_length);
         // check whether the trial iterate is accepted
         bool acceptable_iterate = false;
         if (this->constraint_relaxation_strategy.is_iterate_acceptable(statistics, current_iterate, trial_iterate, direction, step_length)) {
            this->total_number_iterations += this->number_iterations;
            // check termination criteria
            trial_iterate.status = this->check_termination(model, trial_iterate);
            acceptable_iterate = true;
         }
         else if (step_length < this->minimum_step_length) {  // rejected, but small radius
            DEBUG << "The line search step length is smaller than " << this->minimum_step_length << '\n';
            acceptable_iterate = this->terminate_with_small_step(model, direction, trial_iterate);
            if (not acceptable_iterate) {
               //throw std::runtime_error("Line-search strategy reverting to solving the feasibility problem. Not implemented yet.");
               warmstart_information.set_cold_start();
               // revert to solving the feasibility problem
               this->solving_feasibility_problem = true;
               Direction direction_feasibility = this->constraint_relaxation_strategy.solve_feasibility_problem(statistics, current_iterate,
                     direction.primals, warmstart_information);
               BacktrackingLineSearch::check_unboundedness(direction_feasibility);
               Iterate trial_iterate_feasibility = this->backtrack_along_direction(statistics, model, current_iterate, direction_feasibility,
                     warmstart_information);
               this->total_number_iterations += this->number_iterations;
               return trial_iterate_feasibility;
            }
         }

         if (acceptable_iterate) {
            this->set_statistics(statistics, direction, step_length);
            return trial_iterate;
         }
         else {
            step_length = this->decrease_step_length(step_length);
         }
      }
      catch (const EvaluationError& e) {
         WARNING << YELLOW << e.what() << RESET;
         step_length = this->decrease_step_length(step_length);
      }
   }
}

Iterate BacktrackingLineSearch::assemble_trial_iterate(const Model& model, Iterate& current_iterate, const Direction& direction, double step_length) {
   Iterate trial_iterate = GlobalizationMechanism::assemble_trial_iterate(current_iterate, direction, step_length, direction.bound_dual_step_length);

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

bool BacktrackingLineSearch::termination(double primal_dual_step_length) const {
   return (primal_dual_step_length < this->minimum_step_length);
}

void BacktrackingLineSearch::check_unboundedness(const Direction& direction) {
   if (direction.status == SubproblemStatus::UNBOUNDED_PROBLEM) {
      throw std::runtime_error("The subproblem is unbounded, this should not happen. If the subproblem has curvature, use regularization. If not, "
                               "use a trust-region method.\n");
   }
}

void BacktrackingLineSearch::set_statistics(Statistics& statistics, const Direction& direction, double primal_dual_step_length) const {
   statistics.add_statistic("LS iters", this->total_number_iterations);
   statistics.add_statistic("LS step length", primal_dual_step_length);
   statistics.add_statistic("step norm", primal_dual_step_length * direction.norm);
}

void BacktrackingLineSearch::print_iteration(double primal_dual_step_length) {
   DEBUG << "\tLINE SEARCH iteration " << this->number_iterations << ", step_length " << primal_dual_step_length << '\n';
}
