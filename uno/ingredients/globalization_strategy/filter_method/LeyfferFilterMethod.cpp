// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "LeyfferFilterMethod.hpp"

LeyfferFilterMethod::LeyfferFilterMethod(Statistics& /*statistics*/, bool accept_when_switching_violated, const Options& options):
      FilterMethod(options),
      accept_when_switching_violated(accept_when_switching_violated) {
}

/* check acceptability of step(s) (filter & sufficient reduction)
 * filter methods enforce an *unconstrained* sufficient decrease condition
 * precondition: feasible step
 * */
bool LeyfferFilterMethod::is_iterate_acceptable(Statistics& /*statistics*/, const Iterate& /*trial_iterate*/,
      const ProgressMeasures& current_progress_measures, const ProgressMeasures& trial_progress_measures, const ProgressMeasures& predicted_reduction,
      double /*objective_multiplier*/) {
   const double current_optimality_measure = current_progress_measures.optimality(1.) + current_progress_measures.auxiliary_terms;
   const double trial_optimality_measure = trial_progress_measures.optimality(1.) + trial_progress_measures.auxiliary_terms;
   // unconstrained predicted reduction:
   // - ignore the predicted infeasibility reduction
   // - scale the scaled optimality measure with 1
   const double unconstrained_predicted_reduction = predicted_reduction.optimality(1.) + predicted_reduction.auxiliary_terms;
   DEBUG << "Current: η = " << current_progress_measures.infeasibility << ",\t ω = " << current_optimality_measure << '\n';
   DEBUG << "Trial:   η = " << trial_progress_measures.infeasibility << ",\t ω = " << trial_optimality_measure << '\n';
   DEBUG << "Unconstrained predicted reduction: " << predicted_reduction.optimality(1.) << " + " << predicted_reduction.auxiliary_terms <<
         " = " <<  unconstrained_predicted_reduction << '\n';
   DEBUG << *this->filter << '\n';

   bool accept = false;
   // check acceptance
   const bool filter_acceptable = this->filter->acceptable(trial_progress_measures.infeasibility, trial_optimality_measure);
   if (filter_acceptable) {
      DEBUG << "Filter acceptable\n";

      // check acceptance wrt current point
      const bool improves_current_iterate = this->filter->acceptable_wrt_current_iterate(current_progress_measures.infeasibility,
            current_optimality_measure, trial_progress_measures.infeasibility, trial_optimality_measure);
      if (improves_current_iterate) {
         DEBUG << "Acceptable with respect to current point\n";
         const double actual_reduction = this->filter->compute_actual_reduction(current_optimality_measure, current_progress_measures.infeasibility,
               trial_optimality_measure);
         DEBUG << "Actual reduction: " << actual_reduction << '\n';

         // switching condition: the unconstrained predicted reduction is sufficiently positive
         if (this->switching_condition(unconstrained_predicted_reduction, current_progress_measures.infeasibility, this->parameters.delta)) {
            // unconstrained Armijo sufficient decrease condition (predicted reduction should be positive)
            if (this->armijo_sufficient_decrease(unconstrained_predicted_reduction, actual_reduction)) {
               DEBUG << "Trial iterate was accepted by satisfying the Armijo condition\n";
               accept = true;
            }
            else { // switching condition holds, but not Armijo condition
               DEBUG << "Trial iterate was rejected by violating the Armijo condition\n";
            }
         }
         else if (this->accept_when_switching_violated) { // switching condition violated: predicted reduction is not promising
            DEBUG << "Trial iterate was accepted by violating the switching condition\n";
            accept = true;
            DEBUG << "Current iterate was added to the filter\n";
            this->filter->add(current_progress_measures.infeasibility, current_optimality_measure);
         }
         else {
            DEBUG << "Trial iterate was rejected by violating the switching condition\n";
         }
      }
      else {
         DEBUG << "Not acceptable with respect to current point\n";
      }
   }
   else {
      DEBUG << "Not filter acceptable\n";
   }
   DEBUG << '\n';
   return accept;
}