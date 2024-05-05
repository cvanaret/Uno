// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "l1MeritFunction.hpp"

l1MeritFunction::l1MeritFunction(const Options& options): GlobalizationStrategy(options) {
}

void l1MeritFunction::initialize(Statistics& statistics, const Iterate& /*initial_iterate*/, const Options& options) {
   statistics.add_column("penalty param.", Statistics::double_width, options.get_int("statistics_penalty_parameter_column_order"));
}

void l1MeritFunction::reset() {
}

void l1MeritFunction::register_current_progress(const ProgressMeasures& /*current_progress*/) {
}

double l1MeritFunction::get_infeasibility_upper_bound() const {
   // no upper bound monitored
   return INF<double>;
}

void l1MeritFunction::set_infeasibility_upper_bound(double /*new_upper_bound*/, double /*current_infeasibility*/, double /*trial_infeasibility*/) {
   // do nothing
}

bool l1MeritFunction::is_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
      const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction, double objective_multiplier) {
   // predicted reduction with all contributions. This quantity should be positive (= negative directional derivative)
   double constrained_predicted_reduction = l1MeritFunction::constrained_merit_function(predicted_reduction, objective_multiplier);
   DEBUG << "Constrained predicted reduction: " << constrained_predicted_reduction << '\n';
   if (constrained_predicted_reduction <= 0.) {
      WARNING << YELLOW << "The direction is not a descent direction for the merit function. You should decrease the penalty parameter.\n" << RESET;
   }
   // compute current exact penalty
   const double current_merit_value = l1MeritFunction::constrained_merit_function(current_progress, objective_multiplier);
   const double trial_merit_value = l1MeritFunction::constrained_merit_function(trial_progress, objective_multiplier);
   const double actual_reduction = this->compute_merit_actual_reduction(current_merit_value, trial_merit_value);
   DEBUG << "Current merit: " << current_merit_value << '\n';
   DEBUG << "Trial merit:   " << trial_merit_value << '\n';
   DEBUG << "Actual reduction: " << current_merit_value << " - " << trial_merit_value << " = " << actual_reduction << '\n';
   statistics.set("penalty param.", objective_multiplier);

   // Armijo sufficient decrease condition
   const bool accept = this->armijo_sufficient_decrease(constrained_predicted_reduction, actual_reduction);
   if (accept) {
      DEBUG << "Trial iterate was accepted by satisfying Armijo condition\n";
      this->smallest_known_infeasibility = std::min(this->smallest_known_infeasibility, trial_progress.infeasibility);
      statistics.set("status", "accepted (Armijo)");
   }
   else {
      statistics.set("status", "rejected (Armijo)");
   }
   return accept;
}

bool l1MeritFunction::is_feasibility_iterate_acceptable(const ProgressMeasures& /*current_progress*/, const ProgressMeasures& trial_progress) const {
   // if the trial infeasibility improves upon the best known infeasibility
   // TODO put constant in option file
   return (trial_progress.infeasibility <= 0.9*this->smallest_known_infeasibility);
}

double l1MeritFunction::constrained_merit_function(const ProgressMeasures& progress, double objective_multiplier) {
   return progress.objective(objective_multiplier) + progress.auxiliary + progress.infeasibility;
}

double l1MeritFunction::compute_merit_actual_reduction(double current_merit_value, double trial_merit_value) const {
   double actual_reduction = current_merit_value - trial_merit_value;
   if (this->protect_actual_reduction_against_roundoff) {
      static double machine_epsilon = std::numeric_limits<double>::epsilon();
      actual_reduction += 10. * machine_epsilon * std::abs(current_merit_value);
   }
   return actual_reduction;
}