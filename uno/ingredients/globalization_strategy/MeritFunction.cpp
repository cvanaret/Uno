// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "MeritFunction.hpp"

MeritFunction::MeritFunction(const Options& options):
      GlobalizationStrategy(options) {
}

void MeritFunction::initialize(const Iterate& /*first_iterate*/) {
}

void MeritFunction::reset() {
}

void MeritFunction::register_current_progress(const ProgressMeasures& /*current_progress*/) {
}

bool MeritFunction::is_acceptable(const ProgressMeasures& current_progress, const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction) {
   // include the predicted optimality and infeasibility reductions
   const double constrained_predicted_reduction = predicted_reduction.scaled_optimality + predicted_reduction.infeasibility;
   DEBUG << "Constrained predicted reduction: " << constrained_predicted_reduction << '\n';
   // compute current exact penalty
   const double current_exact_merit = current_progress.scaled_optimality + current_progress.unscaled_optimality + current_progress.infeasibility;
   const double trial_exact_merit = trial_progress.scaled_optimality + trial_progress.unscaled_optimality + trial_progress.infeasibility;
   const double actual_reduction = current_exact_merit - trial_exact_merit;
   DEBUG << "Current merit: " << current_progress.scaled_optimality << " + " << current_progress.unscaled_optimality << " + " <<
      current_progress.infeasibility << " = " << current_exact_merit << '\n';
   DEBUG << "Trial merit:   " << trial_progress.scaled_optimality << " + " << trial_progress.unscaled_optimality << " + " <<
      trial_progress.infeasibility << " = " << trial_exact_merit << '\n';
   DEBUG << "Actual reduction: " << current_exact_merit << " - " << trial_exact_merit << " = " << actual_reduction << '\n';

   GlobalizationStrategy::check_finiteness(current_progress);
   GlobalizationStrategy::check_finiteness(trial_progress);

   // Armijo sufficient decrease condition
   const bool accept = this->armijo_sufficient_decrease(constrained_predicted_reduction, actual_reduction);
   if (accept) {
      DEBUG << "Trial iterate was accepted by satisfying Armijo condition\n";
   }
   return accept;
}