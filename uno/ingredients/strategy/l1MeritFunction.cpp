#include "l1MeritFunction.hpp"

/*
 * Infeasibility detection and SQP methods for nonlinear optimization 
 * http://epubs.siam.org/doi/pdf/10.1137/080738222
 */

l1MeritFunction::l1MeritFunction(const Options& options) :
      GlobalizationStrategy(), decrease_fraction(stod(options.at("armijo_decrease_fraction"))) {
}

void l1MeritFunction::initialize(Statistics& statistics, const Iterate& /*first_iterate*/) {
   statistics.add_column("penalty param.", Statistics::double_width, 4);
}

void l1MeritFunction::reset() {
}

void l1MeritFunction::notify(Iterate& /*current_iterate*/) {
}

bool l1MeritFunction::check_acceptance(Statistics& statistics, const ProgressMeasures& current_progress, const ProgressMeasures& trial_progress,
      double objective_multiplier, double predicted_reduction) {
   // compute current exact l1 penalty: rho f + ||c||
   const double current_exact_l1_merit = objective_multiplier * current_progress.objective + current_progress.infeasibility;
   const double trial_exact_l1_merit = objective_multiplier * trial_progress.objective + trial_progress.infeasibility;
   DEBUG << "Current l1 merit: " << current_exact_l1_merit << "\n";
   DEBUG << "Trial l1 merit: " << trial_exact_l1_merit << "\n";
   const double actual_reduction = current_exact_l1_merit - trial_exact_l1_merit;
   DEBUG << "Predicted reduction: " << predicted_reduction << "\n";
   DEBUG << "Actual reduction: " << actual_reduction << "\n";

   bool accept = false;
   // Armijo sufficient decrease condition
   if (actual_reduction >= this->decrease_fraction * predicted_reduction) {
      accept = true;
   }

   if (accept) {
      statistics.add_statistic("penalty param.", objective_multiplier);
      DEBUG << "Trial iterate was accepted by satisfying Armijo condition\n";
   }
   return accept;
}