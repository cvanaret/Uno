// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project root for details.

#include "MeritFunction.hpp"

MeritFunction::MeritFunction(const Options& options):
      GlobalizationStrategy(options) {
}

void MeritFunction::initialize(Statistics& statistics, const Iterate& /*first_iterate*/) {
   statistics.add_column("penalty param.", Statistics::double_width, 4);
}

void MeritFunction::reset() {
}

void MeritFunction::notify(Iterate& /*current_iterate*/) {
}

bool MeritFunction::is_acceptable(Statistics& statistics, const ProgressMeasures& current_progress, const ProgressMeasures& trial_progress,
      double objective_multiplier, double predicted_reduction) {
   GlobalizationStrategy::check_finiteness(current_progress);
   GlobalizationStrategy::check_finiteness(trial_progress);

   DEBUG << "Predicted reduction: " << predicted_reduction << '\n';
   // compute current exact penalty: rho f + ||c||
   const double current_exact_merit = objective_multiplier * current_progress.optimality + current_progress.infeasibility;
   const double trial_exact_merit = objective_multiplier * trial_progress.optimality + trial_progress.infeasibility;
   const double actual_reduction = current_exact_merit - trial_exact_merit;
   DEBUG << "Current merit: " << objective_multiplier << "*" << current_progress.optimality << " + " << current_progress.infeasibility << '\n';
   DEBUG << "Trial merit:   " << objective_multiplier << "*" << trial_progress.optimality << " + " << trial_progress.infeasibility << '\n';
   DEBUG << "Actual reduction: " << current_exact_merit << " - " << trial_exact_merit << " = " << actual_reduction << '\n';

   // Armijo sufficient decrease condition
   const bool accept = this->armijo_sufficient_decrease(predicted_reduction, actual_reduction);
   if (accept) {
      statistics.add_statistic("penalty param.", objective_multiplier);
      DEBUG << "Trial iterate was accepted by satisfying Armijo condition\n";
   }
   return accept;
}