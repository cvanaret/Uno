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
bool WaechterFilterMethod::is_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
      const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction, double objective_multiplier) {
   // in filter methods, we construct an unconstrained measure by ignoring infeasibility and scaling the objective measure by 1
   const double current_objective_measure = current_progress.objective(1.) + current_progress.auxiliary;
   const double trial_objective_measure = trial_progress.objective(1.) + trial_progress.auxiliary;
   const double unconstrained_predicted_reduction = predicted_reduction.objective(1.) + predicted_reduction.auxiliary;
   DEBUG << "Current (infeas., objective+auxiliary) = (" << current_progress.infeasibility << ", " << current_objective_measure << ")\n";
   DEBUG << "Trial   (infeas., objective+auxiliary) = (" << trial_progress.infeasibility << ", " << trial_objective_measure << ")\n";
   DEBUG << "Unconstrained predicted reduction: " << unconstrained_predicted_reduction << '\n';
   DEBUG << "Current filter:\n" << *this->filter;

   const bool solving_feasibility_problem = (objective_multiplier == 0.);
   std::string scenario;
   bool accept = false;
   // solving the feasibility problem = working on infeasibility only (no filter acceptability test)
   if (solving_feasibility_problem) {
      if (this->armijo_sufficient_decrease(predicted_reduction.infeasibility, current_progress.infeasibility - trial_progress.infeasibility)) {
         DEBUG << "Trial iterate (h-type) was accepted by satisfying the Armijo condition\n";
         accept = true;
         if (this->filter->acceptable(trial_progress.infeasibility, trial_objective_measure)) {
            this->filter->add(current_progress.infeasibility, current_objective_measure);
            DEBUG << "Current iterate was added to the filter\n";
         }
      }
      scenario = "h-type Armijo";
   }
   else if (this->filter->acceptable(trial_progress.infeasibility, trial_objective_measure)) {
      DEBUG << "Filter acceptable\n";
      // compute actual reduction
      const double actual_reduction = this->compute_actual_objective_reduction(current_objective_measure, current_progress.infeasibility,
            trial_objective_measure);
      DEBUG << "Actual reduction: " << actual_reduction << '\n';

      // TODO put this coefficient in the option file
      const bool small_infeasibility = current_progress.infeasibility <= 1e-4 * std::max(1., this->initial_infeasibility);
      const bool switching = (0. < unconstrained_predicted_reduction) && this->switching_condition(unconstrained_predicted_reduction,
            current_progress.infeasibility, this->parameters.delta);
      const bool sufficient_decrease = this->armijo_sufficient_decrease(unconstrained_predicted_reduction, actual_reduction);

      // switching condition: the unconstrained predicted reduction is sufficiently positive
      if (small_infeasibility && switching) {
         DEBUG << "Switching condition satisfied\n";
         // unconstrained Armijo sufficient decrease condition: predicted reduction should be positive (f-type)
         if (sufficient_decrease) {
            DEBUG << "Trial iterate (f-type) was accepted by satisfying Armijo condition\n";
            accept = true;
         }
         else {
            DEBUG << "Armijo condition not satisfied\n";
         }
         scenario = "Armijo";
      }
      else {
         DEBUG << "Switching condition violated\n";
         if (this->filter->acceptable_wrt_current_iterate(current_progress.infeasibility, current_objective_measure,
               trial_progress.infeasibility, trial_objective_measure)) {
            DEBUG << "Trial iterate (h-type) acceptable with respect to current point\n";
            accept = true;
         }
         else {
            DEBUG << "Trial iterate (h-type) not acceptable with respect to current point\n";
         }
         scenario = "current point";
      }
      // possibly augment the filter
      if (accept && (not switching || not sufficient_decrease)) {
         DEBUG << "Adding current iterate to the filter\n";
         this->filter->add(current_progress.infeasibility, current_objective_measure);
      }
   }
   else {
      DEBUG << "Trial iterate not filter acceptable\n";
      statistics.set("status", "rejected (filter)");
   }
   statistics.set("status", std::string(accept ? "accepted" : "rejected") + " (" + scenario + ")");
   DEBUG << '\n';
   return accept;
}

bool WaechterFilterMethod::is_feasibility_iterate_acceptable(const ProgressMeasures& current_progress, const ProgressMeasures& trial_progress) const {
   // TODO put constant in the option file
   // TODO current_progress.infeasibility should be replaced with the infeasibility of the first feasibility restoration iterate
   return trial_progress.infeasibility <= 0.9*current_progress.infeasibility &&
      this->filter->acceptable(trial_progress.infeasibility, trial_progress.objective(1.));
}