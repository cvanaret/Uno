// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "WaechterFilterMethod.hpp"
#include "../ProgressMeasures.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Statistics.hpp"
#include "tools/Logger.hpp"

WaechterFilterMethod::WaechterFilterMethod(const Options& options):
      FilterMethod(options) {
}

void WaechterFilterMethod::initialize(Statistics& statistics, const Iterate& initial_iterate, const Options& options) {
   this->initial_infeasibility = initial_iterate.progress.infeasibility;
   FilterMethod::initialize(statistics, initial_iterate, options);
}

/* check acceptability of step(s) (filter & sufficient reduction)
 * filter methods enforce an *unconstrained* sufficient decrease condition
 * precondition: feasible step
 * */
bool WaechterFilterMethod::is_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
      const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction, double objective_multiplier) {
   const bool solving_feasibility_problem = (objective_multiplier == 0.);
   std::string scenario;
   bool accept = false;
   // solving the feasibility problem = working on infeasibility only (no filter acceptability test)
   if (solving_feasibility_problem) {
      DEBUG << "Current infeasibility = " << current_progress.infeasibility << '\n';
      DEBUG << "Trial   infeasibility = " << trial_progress.infeasibility << '\n';
      DEBUG << "Predicted reduction = " << predicted_reduction.infeasibility << '\n';
      if (this->armijo_sufficient_decrease(predicted_reduction.infeasibility, current_progress.infeasibility - trial_progress.infeasibility)) {
         DEBUG << "Trial iterate (h-type) was accepted by satisfying the Armijo condition\n";
         accept = true;
      }
      else {
         DEBUG << "Trial iterate (h-type) was rejected by violating the Armijo condition\n";
      }
      scenario = "h-type Armijo";
      Iterate::number_eval_objective--;
   }
   else {
      // in filter methods, we construct an unconstrained measure by ignoring infeasibility and scaling the objective measure by 1
      const double current_merit = FilterMethod::unconstrained_merit_function(current_progress);
      const double trial_merit = FilterMethod::unconstrained_merit_function(trial_progress);
      const double merit_predicted_reduction = FilterMethod::unconstrained_merit_function(predicted_reduction);
      DEBUG << "Current (infeasibility, objective + auxiliary) = (" << current_progress.infeasibility << ", " << current_merit << ")\n";
      DEBUG << "Trial   (infeasibility, objective + auxiliary) = (" << trial_progress.infeasibility << ", " << trial_merit << ")\n";
      DEBUG << "Current filter:\n" << *this->filter;
      DEBUG << "Unconstrained predicted reduction = " << merit_predicted_reduction << '\n';

      if (this->filter->acceptable(trial_progress.infeasibility, trial_merit)) {
         const double merit_actual_reduction = this->compute_actual_objective_reduction(current_merit, current_progress.infeasibility, trial_merit);
         DEBUG << "Unconstrained actual reduction = " << merit_actual_reduction << '\n';

         // TODO put this coefficient in the option file
         const bool small_infeasibility = current_progress.infeasibility <= 1e-4 * std::max(1., this->initial_infeasibility);
         const bool switching = (0. < merit_predicted_reduction) && this->switching_condition(merit_predicted_reduction, current_progress.infeasibility,
               this->parameters.delta);
         const bool sufficient_decrease = this->armijo_sufficient_decrease(merit_predicted_reduction, merit_actual_reduction);

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
            if (this->filter->acceptable_wrt_current_iterate(current_progress.infeasibility, current_merit, trial_progress.infeasibility, trial_merit)) {
               DEBUG << "Trial iterate (h-type) acceptable with respect to current point\n";
               accept = true;
            }
            else {
               DEBUG << "Trial iterate (h-type) not acceptable with respect to current point\n";
            }
            scenario = "h-type";
         }
         // possibly augment the filter
         if (accept && (not switching || not sufficient_decrease)) {
            DEBUG << "Adding current iterate to the filter\n";
            this->filter->add(current_progress.infeasibility, current_merit);
         }
      }
      else {
         DEBUG << "Trial iterate not filter acceptable\n";
         statistics.set("status", "rejected (filter)");
         scenario = "filter";
      }
   }
   statistics.set("status", std::string(accept ? "accepted" : "rejected") + " (" + scenario + ")");
   return accept;
}

bool WaechterFilterMethod::is_infeasibility_sufficiently_reduced(const ProgressMeasures& current_progress, const ProgressMeasures& trial_progress) const {
   // TODO put constant in the option file
   // TODO current_progress.infeasibility should be replaced with the infeasibility of the first feasibility restoration iterate
   return trial_progress.infeasibility <= 0.9 * current_progress.infeasibility &&
      this->filter->acceptable(trial_progress.infeasibility, FilterMethod::unconstrained_merit_function(trial_progress));
}