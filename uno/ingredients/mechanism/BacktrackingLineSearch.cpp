#include <cmath>
#include <cassert>
#include "BacktrackingLineSearch.hpp"
#include "tools/Logger.hpp"

BacktrackingLineSearch::BacktrackingLineSearch(ConstraintRelaxationStrategy& constraint_relaxation_strategy, int max_iterations, double backtracking_ratio):
   GlobalizationMechanism(constraint_relaxation_strategy, max_iterations), backtracking_ratio(backtracking_ratio) {
}

void BacktrackingLineSearch::initialize(Statistics& statistics, const Problem& problem, Iterate& first_iterate) {
   statistics.add_column("SOC", Statistics::char_width, 9);
   statistics.add_column("LS step length", Statistics::double_width, 30);

   // generate the initial point
   this->relaxation_strategy.initialize(statistics, problem, first_iterate);
}

std::tuple<Iterate, double> BacktrackingLineSearch::compute_acceptable_iterate(Statistics& statistics, const Problem& problem, Iterate&
current_iterate) {
   /* compute the direction */
   this->relaxation_strategy.generate_subproblem(problem, current_iterate, INFINITY);
   Direction direction = this->relaxation_strategy.compute_feasible_direction(statistics, problem, current_iterate);

   /* step length follows the following sequence: 1, ratio, ratio^2, ratio^3, ... */
   this->step_length = 1.;
   this->number_iterations = 0;
   bool failure = false;
   while (!failure) {
      while (!this->termination_()) {
         assert(0 < this->step_length && this->step_length <= 1);
         this->number_iterations++;
         this->print_iteration_();

         // assemble the trial iterate
         Iterate trial_iterate = this->assemble_trial_iterate(current_iterate, direction, this->step_length);

         try {
            // check whether the trial step is accepted
            if (this->relaxation_strategy.is_acceptable(statistics, problem, current_iterate, trial_iterate, direction, this->step_length)) {
               this->add_statistics(statistics, direction);

               // let the subproblem know the accepted iterate
               this->relaxation_strategy.register_accepted_iterate(trial_iterate);
               return std::make_tuple(std::move(trial_iterate), direction.norm);
            }
            else if (this->number_iterations == 1 && trial_iterate.progress.infeasibility >= current_iterate.progress.infeasibility) { // reject
               // compute a (temporary) SOC direction
               Direction direction_soc = this->relaxation_strategy.compute_second_order_correction(problem, trial_iterate);

               // assemble the (temporary) SOC trial iterate
               Iterate trial_iterate_soc = this->assemble_trial_iterate(current_iterate, direction_soc, this->step_length);

               if (this->relaxation_strategy.is_acceptable(statistics, problem, current_iterate, trial_iterate_soc, direction_soc,
                     this->step_length)) {
                  this->add_statistics(statistics, direction_soc);
                  statistics.add_statistic("SOC", "x");

                  // let the subproblem know the accepted iterate
                  this->relaxation_strategy.register_accepted_iterate(trial_iterate_soc);
                  trial_iterate_soc.multipliers.lower_bounds = trial_iterate.multipliers.lower_bounds;
                  trial_iterate_soc.multipliers.upper_bounds = trial_iterate.multipliers.upper_bounds;
                  return std::make_tuple(std::move(trial_iterate_soc), direction_soc.norm);
               }
               else {
                  /* decrease the step length */
                  DEBUG << "SOC step discarded\n\n";
                  statistics.add_statistic("SOC", "-");
                  this->decrease_step_length();
               }
            }
            else {
               this->decrease_step_length();
            }
         }
         catch (const NumericalError& e) {
            GlobalizationMechanism::print_warning(e.what());
            this->decrease_step_length();
         }
      }
      // if step length is too small, run restoration phase
      if (this->step_length < this->min_step_length && 0. < direction.multipliers.objective) {
         assert(false);
         //if (0. < current_iterate.progress.feasibility && !direction.is_relaxed) {
         // reset the line search with the restoration solution
         direction = this->relaxation_strategy.solve_feasibility_problem(statistics, problem, current_iterate, direction);
         this->step_length = 1.;
         this->number_iterations = 0;
      }
      else {
         failure = true;
      }
   }
   throw std::runtime_error("Line search: maximum number of iterations reached");
}

void BacktrackingLineSearch::add_statistics(Statistics& statistics, const Direction& direction) {
   statistics.add_statistic("minor", this->number_iterations);
   statistics.add_statistic("LS step length", this->step_length);
   statistics.add_statistic("step norm", this->step_length * direction.norm);
}

void BacktrackingLineSearch::decrease_step_length() {
   this->step_length *= this->backtracking_ratio;
}

bool BacktrackingLineSearch::termination_() {
   return (this->max_iterations < this->number_iterations);
}

void BacktrackingLineSearch::print_iteration_() {
   DEBUG << "\n\tLINE SEARCH iteration " << this->number_iterations << ", step_length " << this->step_length << "\n";
}