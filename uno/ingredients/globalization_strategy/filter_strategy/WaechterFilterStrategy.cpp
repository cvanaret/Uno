// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "WaechterFilterStrategy.hpp"

WaechterFilterStrategy::WaechterFilterStrategy(Statistics& /*statistics*/, const Options& options): FilterStrategy(options) {
}

void WaechterFilterStrategy::initialize(const Iterate& initial_iterate) {
   this->initial_infeasibility = initial_iterate.residuals.infeasibility;
   FilterStrategy::initialize(initial_iterate);
}

/* check acceptability of step(s) (filter & sufficient reduction)
 * filter methods enforce an *unconstrained* sufficient decrease condition
 * precondition: feasible step
 * */
bool WaechterFilterStrategy::is_iterate_acceptable(Statistics& /*statistics*/, const Iterate& /*trial_iterate*/,
      const ProgressMeasures& current_progress_measures, const ProgressMeasures& trial_progress_measures, const ProgressMeasures& predicted_reduction,
      double /*objective_multiplier*/) {
   const double current_optimality_measure = current_progress_measures.optimality(1.) + current_progress_measures.auxiliary_terms;
   const double trial_optimality_measure = trial_progress_measures.optimality(1.) + trial_progress_measures.auxiliary_terms;
   // unconstrained predicted reduction:
   // - ignore the predicted infeasibility reduction
   // - scale the scaled optimality measure with 1
   const double unconstrained_predicted_reduction = predicted_reduction.optimality(1.) + predicted_reduction.auxiliary_terms;
   DEBUG << "Current: η = " << current_progress_measures.infeasibility << ", ω = " << current_progress_measures.optimality(1.) << " + " <<
      current_progress_measures.auxiliary_terms << '\n';
   DEBUG << "Trial:   η = " << trial_progress_measures.infeasibility << ", ω = " << trial_progress_measures.optimality(1.) << " + " <<
      trial_progress_measures.auxiliary_terms << '\n';
   DEBUG << "Unconstrained predicted reduction: " << predicted_reduction.optimality(1.) << " + " << predicted_reduction.auxiliary_terms <<
         " = " <<  unconstrained_predicted_reduction << '\n';

   GlobalizationStrategy::check_finiteness(current_progress_measures, 1.);
   GlobalizationStrategy::check_finiteness(trial_progress_measures, 1.);
   DEBUG << *this->filter;

   bool accept = false;
   // check acceptance
   const bool filter_acceptable = this->filter->acceptable(trial_progress_measures.infeasibility, trial_optimality_measure);
   if (filter_acceptable) {
      DEBUG << "Filter acceptable\n";
      // compute actual reduction (and protect against roundoff errors)
      static double machine_epsilon = std::numeric_limits<double>::epsilon();
      const double actual_reduction = this->filter->compute_actual_reduction(current_optimality_measure, current_progress_measures.infeasibility,
            trial_optimality_measure) + 10. * machine_epsilon * std::abs(current_optimality_measure);
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
         }
         else {
            DEBUG << "Armijo condition not satisfied\n";
         }
      }
      else {
         DEBUG << "Switching condition violated\n";
         if (this->filter->acceptable_wrt_current_iterate(current_progress_measures.infeasibility, current_optimality_measure,
               trial_progress_measures.infeasibility, trial_optimality_measure)) {
            DEBUG << "Acceptable wrt current point\n";
            accept = true;
         }
         else {
            DEBUG << "Not acceptable wrt current point\n";
         }
      }
      // possibly augment the filter
      if (accept && (not switching || not sufficient_decrease)) {
         DEBUG << "Adding current iterate to the filter\n";
         this->filter->add(current_progress_measures.infeasibility, current_optimality_measure);
      }
   }
   else {
      DEBUG << "Not filter acceptable\n";
   }
   DEBUG << '\n';
   return accept;
}
