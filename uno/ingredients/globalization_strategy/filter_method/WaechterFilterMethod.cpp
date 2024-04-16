// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "WaechterFilterMethod.hpp"

WaechterFilterMethod::WaechterFilterMethod(const Options& options):
      FilterMethod(options) {
}

void WaechterFilterMethod::initialize(Statistics& statistics, const Iterate& initial_iterate, const Options& options) {
   this->initial_infeasibility = initial_iterate.residuals.infeasibility;
   FilterMethod::initialize(statistics, initial_iterate, options);
}

/* check acceptability of step(s) (filter & sufficient reduction)
 * filter methods enforce an *unconstrained* sufficient decrease condition
 * precondition: feasible step
 * */
bool WaechterFilterMethod::is_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress_measures,
      const ProgressMeasures& trial_progress_measures, const ProgressMeasures& predicted_reduction, double /*objective_multiplier*/) {
   const double current_objective_measure = current_progress_measures.objective(1.) + current_progress_measures.auxiliary;
   const double trial_objective_measure = trial_progress_measures.objective(1.) + trial_progress_measures.auxiliary;

   // unconstrained predicted reduction:
   // - ignore the predicted infeasibility reduction
   // - scale the objective measure with 1
   const double unconstrained_predicted_reduction = predicted_reduction.objective(1.) + predicted_reduction.auxiliary;
   DEBUG << "Unconstrained predicted reduction: " << unconstrained_predicted_reduction << '\n';
   DEBUG << "Current filter:\n" << *this->filter;

   bool accept = false;
   // check acceptance
   const bool filter_acceptable = this->filter->acceptable(trial_progress_measures.infeasibility, trial_objective_measure);
   if (filter_acceptable) {
      DEBUG << "Filter acceptable\n";
      DEBUG << "Current (infeas., objective+auxiliary) = (" << current_progress_measures.infeasibility << ", " << current_objective_measure << ")\n";
      DEBUG << "Trial   (infeas., objective+auxiliary) = (" << trial_progress_measures.infeasibility << ", " << trial_objective_measure << ")\n";

      // compute actual reduction (and protect against roundoff errors)
      // TODO put constants in the option file
      static double machine_epsilon = std::numeric_limits<double>::epsilon();
      const double actual_reduction = this->filter->compute_actual_reduction(current_objective_measure, current_progress_measures.infeasibility,
            trial_objective_measure) + 10. * machine_epsilon * std::abs(current_objective_measure);
      DEBUG << "Actual reduction: " << actual_reduction << '\n';

      // TODO put this coefficient in the option file
      const bool small_infeasibility = current_progress_measures.infeasibility <= 1e-4*std::max(1., this->initial_infeasibility);
      const bool switching = (0. < unconstrained_predicted_reduction) && this->switching_condition(unconstrained_predicted_reduction,
            current_progress_measures.infeasibility, this->parameters.delta);
      const bool sufficient_decrease = this->armijo_sufficient_decrease(unconstrained_predicted_reduction, actual_reduction);

      // switching condition: the unconstrained predicted reduction is sufficiently positive
      if (small_infeasibility && switching) {
         DEBUG << "Switching condition satisfied\n";
         // unconstrained Armijo sufficient decrease condition (predicted reduction should be positive)
         if (sufficient_decrease) {
            DEBUG << "Trial iterate was accepted by satisfying Armijo condition\n";
            accept = true;
            statistics.set("status", "accepted (Armijo)");
         }
         else {
            DEBUG << "Armijo condition not satisfied\n";
            statistics.set("status", "rejected (Armijo)");
         }
      }
      else {
         DEBUG << "Switching condition violated\n";
         if (this->filter->acceptable_wrt_current_iterate(current_progress_measures.infeasibility, current_objective_measure,
               trial_progress_measures.infeasibility, trial_objective_measure)) {
            DEBUG << "Acceptable with respect to current point\n";
            accept = true;
            statistics.set("status", "accepted (current point)");
         }
         else {
            DEBUG << "Not acceptable with respect to current point\n";
            statistics.set("status", "rejected (current point)");
         }
      }
      // possibly augment the filter
      if (accept && (not switching || not sufficient_decrease)) {
         DEBUG << "Adding current iterate to the filter\n";
         this->filter->add(current_progress_measures.infeasibility, current_objective_measure);
      }
   }
   else {
      DEBUG << "Not filter acceptable\n";
      statistics.set("status", "rejected (filter)");
   }
   DEBUG << '\n';
   return accept;
}
