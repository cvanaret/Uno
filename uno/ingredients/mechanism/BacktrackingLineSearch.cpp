#include <cmath>
#include <cassert>
#include "ingredients/subproblem/AugmentedSystem.hpp"
#include "BacktrackingLineSearch.hpp"
#include "tools/Logger.hpp"

BacktrackingLineSearch::BacktrackingLineSearch(ConstraintRelaxationStrategy& constraint_relaxation_strategy, const Options& options):
      GlobalizationMechanism(constraint_relaxation_strategy),
      backtracking_ratio(std::stod(options.at("LS_backtracking_ratio"))),
      min_step_length(std::stod(options.at("LS_min_step_length"))),
      use_second_order_correction(options.at("use_second_order_correction") == "yes") {
   assert(0 < this->backtracking_ratio && this->backtracking_ratio < 1. && "The LS backtracking ratio should be in (0, 1)");
   assert(0 < this->min_step_length && this->min_step_length < 1. && "The LS minimum step length should be in (0, 1)");
}

void BacktrackingLineSearch::initialize(Statistics& statistics, Iterate& first_iterate) {
   statistics.add_column("SOC", Statistics::char_width, 9);
   statistics.add_column("LS step length", Statistics::double_width, 30);

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
      return this->constraint_relaxation_strategy.solve_feasibility_problem(statistics, current_iterate, std::nullopt);
   }
}

std::tuple<Iterate, double> BacktrackingLineSearch::compute_acceptable_iterate(Statistics& statistics, Iterate& current_iterate) {
   // compute the direction
   this->constraint_relaxation_strategy.set_variable_bounds(current_iterate, std::numeric_limits<double>::infinity());
   Direction direction = this->compute_direction(statistics, current_iterate);
   GlobalizationMechanism::check_unboundedness(direction);
   PredictedReductionModel predicted_reduction_model = this->constraint_relaxation_strategy.generate_predicted_reduction_model(direction);
   this->solving_feasibility_problem = false;

   // step length follows the following sequence: 1, ratio, ratio^2, ratio^3, ...
   this->step_length = 1.;
   this->number_iterations = 0;
   bool failure = false;
   while (!failure) {
      while (!this->termination()) {
         assert(0 < this->step_length && this->step_length <= 1 && "The line-search step length is not in (0, 1]");
         this->number_iterations++;
         this->print_iteration();

         Iterate trial_iterate = GlobalizationMechanism::assemble_trial_iterate(current_iterate, direction, this->step_length);
         try {
            const bool is_acceptable = this->constraint_relaxation_strategy.is_acceptable(statistics, current_iterate, trial_iterate,
                  direction, predicted_reduction_model, this->step_length);
            // check whether the trial step is accepted
            if (is_acceptable) {
               DEBUG << "Trial step accepted\n\n";
               // let the subproblem know the accepted iterate
               this->constraint_relaxation_strategy.register_accepted_iterate(trial_iterate);
               this->add_statistics(statistics, direction);
               return std::make_tuple(std::move(trial_iterate), direction.norm);
            }
            else if (this->use_second_order_correction && this->constraint_relaxation_strategy.soc_strategy() == SOC_UPON_REJECTION &&
                     this->number_iterations == 1 && trial_iterate.nonlinear_progress.infeasibility >= current_iterate.nonlinear_progress.infeasibility &&
                     !this->solving_feasibility_problem) {
               // reject the full step: compute a (temporary) SOC direction
               Direction direction_soc = this->constraint_relaxation_strategy.compute_second_order_correction(trial_iterate);

               // assemble the (temporary) SOC trial iterate
               Iterate trial_iterate_soc = GlobalizationMechanism::assemble_trial_iterate(current_iterate, direction_soc, this->step_length);

               if (this->constraint_relaxation_strategy.is_acceptable(statistics, current_iterate, trial_iterate_soc, direction_soc,
                     predicted_reduction_model, this->step_length)) {
                  DEBUG << "Trial SOC step accepted\n";
                  this->add_statistics(statistics, direction_soc);
                  statistics.add_statistic("SOC", "x");

                  // let the subproblem know the accepted iterate
                  this->constraint_relaxation_strategy.register_accepted_iterate(trial_iterate_soc);
                  trial_iterate_soc.multipliers.lower_bounds = trial_iterate.multipliers.lower_bounds;
                  trial_iterate_soc.multipliers.upper_bounds = trial_iterate.multipliers.upper_bounds;
                  return std::make_tuple(std::move(trial_iterate_soc), direction_soc.norm);
               }
               else {
                  DEBUG << "SOC step discarded\n\n";
                  statistics.add_statistic("SOC", "-");
                  this->decrease_step_length();
               }
            }
            else { // trial iterate not acceptable
               this->decrease_step_length();
            }
         }
         catch (const NumericalError& e) {
            GlobalizationMechanism::print_warning(e.what());
            this->decrease_step_length();
         }
      }
      // if step length is too small, revert to solving the feasibility problem (if we aren't already solving it)
      if (!this->solving_feasibility_problem && 0. < direction.multipliers.objective) {
         // TODO: test if 0. < current_iterate.progress.feasibility ?
         DEBUG << "The line search failed, switching to feasibility problem\n";
         // reset the line search with the restoration solution
         direction = this->constraint_relaxation_strategy.solve_feasibility_problem(statistics, current_iterate, direction.primals);
         BacktrackingLineSearch::check_unboundedness(direction);
         this->step_length = 1.;
         this->number_iterations = 0;
         this->solving_feasibility_problem = true;
      }
      else {
         WARNING << "The feasibility problem failed to make progress\n";
         failure = true;
      }
   }
   throw std::runtime_error("Line search: maximum number of iterations reached");
}

void BacktrackingLineSearch::decrease_step_length() {
   this->step_length *= this->backtracking_ratio;
}

bool BacktrackingLineSearch::termination() const {
   return (this->step_length < this->min_step_length);
}

void BacktrackingLineSearch::add_statistics(Statistics& statistics, const Direction& direction) {
   statistics.add_statistic("minor", this->number_iterations);
   statistics.add_statistic("LS step length", this->step_length);
   statistics.add_statistic("step norm", this->step_length * direction.norm);
}

void BacktrackingLineSearch::print_iteration() {
   DEBUG << "\tLINE SEARCH iteration " << this->number_iterations << ", step_length " << this->step_length << "\n";
}
