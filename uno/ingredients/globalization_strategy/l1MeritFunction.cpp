// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "l1MeritFunction.hpp"

l1MeritFunction::l1MeritFunction(Statistics& statistics, const Options& options): GlobalizationStrategy(options) {
   statistics.add_column("penalty param.", Statistics::double_width, options.get_int("statistics_penalty_parameter_column_order"));
}

void l1MeritFunction::initialize(const Iterate& /*initial_iterate*/) {
}

void l1MeritFunction::reset() {
}

void l1MeritFunction::register_current_progress(const ProgressMeasures& /*current_progress*/) {
}

bool l1MeritFunction::is_iterate_acceptable(Statistics& statistics, const Iterate& /*trial_iterate*/, const ProgressMeasures& current_progress,
      const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction, double objective_multiplier) {
   // predicted reduction with all contributions. This quantity should be positive (= negative directional derivative)
   double constrained_predicted_reduction = predicted_reduction.optimality(objective_multiplier) + predicted_reduction.auxiliary_terms +
         predicted_reduction.infeasibility;
   DEBUG << "Constrained predicted reduction: " << constrained_predicted_reduction << '\n';
   if (constrained_predicted_reduction <= 0.) {
      WARNING << YELLOW << "The direction is not a descent direction for the merit function. You should decrease the penalty parameter.\n" << RESET;
   }
   // compute current exact penalty
   const double current_exact_merit = current_progress.optimality(objective_multiplier) + current_progress.auxiliary_terms + current_progress.infeasibility;
   const double trial_exact_merit = trial_progress.optimality(objective_multiplier) + trial_progress.auxiliary_terms + trial_progress.infeasibility;
   const double actual_reduction = current_exact_merit - trial_exact_merit;
   DEBUG << "Current merit: " << current_progress.optimality(objective_multiplier) << " + " << current_progress.auxiliary_terms << " + " <<
         current_progress.infeasibility << " = " << current_exact_merit << '\n';
   DEBUG << "Trial merit:   " << trial_progress.optimality(objective_multiplier) << " + " << trial_progress.auxiliary_terms << " + " <<
         trial_progress.infeasibility << " = " << trial_exact_merit << '\n';
   DEBUG << "Actual reduction: " << current_exact_merit << " - " << trial_exact_merit << " = " << actual_reduction << '\n';
   statistics.add_statistic("penalty param.", objective_multiplier);

   // Armijo sufficient decrease condition
   const bool accept = this->armijo_sufficient_decrease(constrained_predicted_reduction, actual_reduction);
   if (accept) {
      DEBUG << "Trial iterate was accepted by satisfying Armijo condition\n";
      this->smallest_known_infeasibility = std::min(this->smallest_known_infeasibility, trial_progress.infeasibility);
   }
   return accept;
}

bool l1MeritFunction::is_infeasibility_acceptable(double infeasibility_measure) const {
   // accept if the infeasibility measure improves upon the smallest known infeasibility
   return (infeasibility_measure < this->smallest_known_infeasibility);
}