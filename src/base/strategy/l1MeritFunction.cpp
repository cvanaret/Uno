#include <cmath>
#include "Uno.hpp"
#include "l1MeritFunction.hpp"
#include "Logger.hpp"

/*
 * Infeasibility detection and SQP methods for nonlinear optimization 
 * http://epubs.siam.org/doi/pdf/10.1137/080738222
 */

l1MeritFunction::l1MeritFunction() : GlobalizationStrategy(), decrease_fraction_(1e-8) {
}

void l1MeritFunction::initialize(Statistics& statistics, const Iterate& /*first_iterate*/) {
   statistics.add_column("penalty param.", Statistics::double_width, 4);
}

void l1MeritFunction::reset() {
}

void l1MeritFunction::notify(Iterate& /*current_iterate*/) {
}

bool l1MeritFunction::check_acceptance(Statistics& statistics, ProgressMeasures& current_progress, ProgressMeasures& trial_progress,
      const Direction& direction, double predicted_reduction) {
   bool accept = false;
   /* compute current exact l1 penalty: rho f + ||c|| */
   double current_exact_l1_penalty = direction.objective_multiplier * current_progress.objective + current_progress.feasibility;
   double trial_exact_l1_penalty = direction.objective_multiplier * trial_progress.objective + trial_progress.feasibility;

   double actual_reduction = current_exact_l1_penalty - trial_exact_l1_penalty;
   DEBUG << "Predicted reduction: " << predicted_reduction << ", actual: " << actual_reduction << "\n\n";
   // Armijo sufficient decrease condition
   if (actual_reduction >= this->decrease_fraction_ * predicted_reduction) {
      accept = true;
   }

   if (accept) {
      statistics.add_statistic("penalty param.", direction.objective_multiplier);
   }
   return accept;
}