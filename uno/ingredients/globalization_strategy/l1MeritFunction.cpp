#include "l1MeritFunction.hpp"

l1MeritFunction::l1MeritFunction(const Options& options):
      GlobalizationStrategy(options) {
}

void l1MeritFunction::initialize(Statistics& statistics, const Iterate& /*first_iterate*/) {
   statistics.add_column("penalty param.", Statistics::double_width, 4);
}

void l1MeritFunction::reset() {
}

void l1MeritFunction::notify(Iterate& /*current_iterate*/) {
}

bool l1MeritFunction::is_acceptable(Statistics& statistics, const ProgressMeasures& current_progress, const ProgressMeasures& trial_progress,
      double objective_multiplier, double predicted_reduction) {
   GlobalizationStrategy::check_finiteness(current_progress);
   GlobalizationStrategy::check_finiteness(trial_progress);

   DEBUG << "Predicted reduction: " << predicted_reduction << '\n';
   // compute current exact l1 penalty: rho f + ||c||
   const double current_exact_l1_merit = objective_multiplier * current_progress.reformulation_objective + current_progress.infeasibility;
   const double trial_exact_l1_merit = objective_multiplier * trial_progress.reformulation_objective + trial_progress.infeasibility;
   const double actual_reduction = current_exact_l1_merit - trial_exact_l1_merit;
   DEBUG << "Current l1 merit: " << objective_multiplier << "*" << current_progress.reformulation_objective << " + " << current_progress.infeasibility << '\n';
   DEBUG << "Trial l1 merit:   " << objective_multiplier << "*" << trial_progress.reformulation_objective << " + " << trial_progress.infeasibility << '\n';
   DEBUG << "Actual reduction: " << current_exact_l1_merit << " - " << trial_exact_l1_merit << " = " << actual_reduction << '\n';

   bool accept = false;
   // Armijo sufficient decrease condition
   if (this->armijo_sufficient_decrease(predicted_reduction, actual_reduction)) {
      accept = true;
   }

   if (accept) {
      statistics.add_statistic("penalty param.", objective_multiplier);
      DEBUG << "Trial iterate was accepted by satisfying Armijo condition\n";
   }
   return accept;
}