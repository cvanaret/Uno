// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "LeyfferFilterMethod.hpp"

LeyfferFilterMethod::LeyfferFilterMethod(bool is_solving_feasibility_problem, const Options& options):
      FilterMethod(options),
      is_solving_feasibility_problem(is_solving_feasibility_problem) {
}

/* check acceptability of step(s) (filter & sufficient reduction)
 * filter methods enforce an *unconstrained* sufficient decrease condition
 * precondition: feasible step
 * */
bool LeyfferFilterMethod::is_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress_measures,
      const ProgressMeasures& trial_progress_measures, const ProgressMeasures& predicted_reduction, double /*objective_multiplier*/) {
   const double current_objective_measure = current_progress_measures.objective(1.) + current_progress_measures.auxiliary;
   const double trial_objective_measure = trial_progress_measures.objective(1.) + trial_progress_measures.auxiliary;

   // unconstrained predicted reduction:
   // - ignore the predicted infeasibility reduction
   // - scale the objective measure with 1
   const double unconstrained_predicted_reduction = predicted_reduction.objective(1.) + predicted_reduction.auxiliary;
   DEBUG << "Unconstrained predicted reduction: " << unconstrained_predicted_reduction << '\n';
   DEBUG << "Current filter:\n" << *this->filter << '\n';

   bool accept = false;
   // check acceptance
   const bool filter_acceptable = this->filter->acceptable(trial_progress_measures.infeasibility, trial_objective_measure);
   if (filter_acceptable) {
      DEBUG << "Filter acceptable\n";
      DEBUG << "Current: (infeas., objective+auxiliary) = (" << current_progress_measures.infeasibility << ", " << current_objective_measure << ")\n";
      DEBUG << "Trial:   (infeas., objective+auxiliary) = (" << trial_progress_measures.infeasibility << ", " << trial_objective_measure << ")\n";

      // check acceptance wrt current point
      const bool improves_current_iterate = this->filter->acceptable_wrt_current_iterate(current_progress_measures.infeasibility,
            current_objective_measure, trial_progress_measures.infeasibility, trial_objective_measure);
      if (improves_current_iterate) {
         DEBUG << "Acceptable with respect to current point\n";
         const double actual_reduction = this->filter->compute_actual_reduction(current_objective_measure, current_progress_measures.infeasibility,
               trial_objective_measure);
         DEBUG << "Actual reduction: " << actual_reduction << '\n';

         // switching condition: the unconstrained predicted reduction is sufficiently positive
         if (this->switching_condition(unconstrained_predicted_reduction, current_progress_measures.infeasibility, this->parameters.delta)) {
            // unconstrained Armijo sufficient decrease condition: predicted reduction should be positive (f-type)
            if (this->armijo_sufficient_decrease(unconstrained_predicted_reduction, actual_reduction)) {
               DEBUG << "Trial iterate (f-type) was accepted by satisfying the Armijo condition\n";
               accept = true;
               statistics.set("status", "accepted (Armijo)");
            }
            else { // switching condition holds, but not Armijo condition
               DEBUG << "Trial iterate was rejected by violating the Armijo condition\n";
               statistics.set("status", "rejected (Armijo)");
            }
         }
         // TODO: switching condition always satisfied in feasibility restoration (pred > 0 and infeasibility = 0). Why is there this test?
         else if (not this->is_solving_feasibility_problem) { // switching condition violated: predicted reduction is not promising (h-type)
            DEBUG << "Trial iterate (h-type) was accepted by violating the switching condition\n";
            accept = true;
            DEBUG << "Current iterate was added to the filter\n";
            this->filter->add(current_progress_measures.infeasibility, current_objective_measure);
            statistics.set("status", "accepted (!switching)");
         }
         else {
            DEBUG << "Trial iterate was rejected by violating the switching condition\n";
            statistics.set("status", "rejected (switching)");
         }
      }
      else {
         DEBUG << "Not acceptable with respect to current point\n";
         statistics.set("status", "rejected (current point)");
      }
   }
   else {
      DEBUG << "Not filter acceptable\n";
      statistics.set("status", "rejected (filter)");
   }
   DEBUG << '\n';
   return accept;
}

bool LeyfferFilterMethod::is_infeasibility_acceptable(const ProgressMeasures& /*current_progress*/, const ProgressMeasures& trial_progress) const {
   // if the trial infeasibility improves upon the best known infeasibility
   return (trial_progress.infeasibility < this->filter->get_smallest_infeasibility());
}