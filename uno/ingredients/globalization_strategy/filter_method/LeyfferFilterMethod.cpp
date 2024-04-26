// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "LeyfferFilterMethod.hpp"

LeyfferFilterMethod::LeyfferFilterMethod(const Options& options): FilterMethod(options) {
}

/* check acceptability of step(s) (filter & sufficient reduction)
 * filter methods enforce an *unconstrained* sufficient decrease condition
 * precondition: feasible step
 * */
bool LeyfferFilterMethod::is_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
      const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction, double objective_multiplier) {
   // in filter methods, we construct an unconstrained measure by ignoring infeasibility and scaling the objective measure by 1
   const double current_objective_measure = current_progress.objective(1.) + current_progress.auxiliary;
   const double trial_objective_measure = trial_progress.objective(1.) + trial_progress.auxiliary;
   const double objective_predicted_reduction = predicted_reduction.objective(1.) + predicted_reduction.auxiliary;
   DEBUG << "Current: (infeas., objective+auxiliary) = (" << current_progress.infeasibility << ", " << current_objective_measure << ")\n";
   DEBUG << "Trial:   (infeas., objective+auxiliary) = (" << trial_progress.infeasibility << ", " << trial_objective_measure << ")\n";
   DEBUG << "Unconstrained predicted reduction: " << objective_predicted_reduction << '\n';
   DEBUG << "Current filter:\n" << *this->filter << '\n';

   const bool solving_feasibility_problem = (objective_multiplier == 0.);
   std::string scenario;
   bool accept = false;
   // solving the feasibility problem = working on infeasibility only (no filter acceptability test)
   if (solving_feasibility_problem) {
      if (this->armijo_sufficient_decrease(predicted_reduction.infeasibility, current_progress.infeasibility - trial_progress.infeasibility)) {
         DEBUG << "Trial iterate (h-type) was accepted by satisfying the Armijo condition\n";
         accept = true;
      }
      else {
         DEBUG << "Trial iterate (h-type) was rejected by violating the Armijo condition\n";
      }
      scenario = "h-type Armijo";
   }
   else if (this->filter->acceptable(trial_progress.infeasibility, trial_objective_measure)) {
      if (this->filter->acceptable_wrt_current_iterate(current_progress.infeasibility, current_objective_measure, trial_progress.infeasibility,
            trial_objective_measure)) {
         // switching condition: check whether the unconstrained predicted reduction is sufficiently positive
         if (this->switching_condition(objective_predicted_reduction, current_progress.infeasibility, this->parameters.delta)) {
            // unconstrained Armijo sufficient decrease condition: predicted reduction should be positive (f-type)
            const double objective_actual_reduction = this->compute_actual_objective_reduction(current_objective_measure, current_progress.infeasibility,
                  trial_objective_measure);
            DEBUG << "Actual reduction: " << objective_actual_reduction << '\n';
            if (this->armijo_sufficient_decrease(objective_predicted_reduction, objective_actual_reduction)) {
               DEBUG << "Trial iterate (f-type) was accepted by satisfying the Armijo condition\n";
               accept = true;
            }
            else { // switching condition holds, but not Armijo condition
               DEBUG << "Trial iterate (f-type) was rejected by violating the Armijo condition\n";
            }
            scenario = "f-type Armijo";
         }
         // switching condition violated: predicted reduction is not promising (h-type)
         else {
            DEBUG << "Trial iterate (h-type) was accepted by violating the switching condition\n";
            accept = true;
            this->filter->add(current_progress.infeasibility, current_objective_measure);
            DEBUG << "Current iterate was added to the filter\n";
            scenario = "!switching";
         }
      }
      else {
         scenario = "current point";
      }
   }
   else {
      scenario = "filter";
   }
   statistics.set("status", std::string(accept ? "accepted" : "rejected") + " (" + scenario + ")");
   DEBUG << '\n';
   return accept;
}

bool LeyfferFilterMethod::is_feasibility_iterate_acceptable(const ProgressMeasures& /*current_progress*/, const ProgressMeasures& trial_progress) const {
   // if the trial infeasibility improves upon the best known infeasibility
   return (trial_progress.infeasibility < this->filter->get_smallest_infeasibility());
}