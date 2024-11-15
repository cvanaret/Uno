// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "FletcherFilterMethod.hpp"
#include "filters/Filter.hpp"
#include "ingredients/globalization_strategies/ProgressMeasures.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Logger.hpp"
#include "options/Options.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   FletcherFilterMethod::FletcherFilterMethod(const Options& options): FilterMethod(options) { }

   FletcherFilterMethod::~FletcherFilterMethod() { }

   bool FletcherFilterMethod::is_regular_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
         const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction) {
      // in filter methods, we construct an unconstrained measure by ignoring infeasibility and scaling the objective measure by 1
      const double current_merit = FilterMethod::unconstrained_merit_function(current_progress);
      const double trial_merit = FilterMethod::unconstrained_merit_function(trial_progress);
      const double merit_predicted_reduction = FilterMethod::unconstrained_merit_function(predicted_reduction);
      DEBUG << "Current: (infeasibility, objective + auxiliary) = (" << current_progress.infeasibility << ", " << current_merit << ")\n";
      DEBUG << "Trial:   (infeasibility, objective + auxiliary) = (" << trial_progress.infeasibility << ", " << trial_merit << ")\n";
      DEBUG << "Current filter:\n" << *this->filter << '\n';
      DEBUG << "Unconstrained predicted reduction = " << merit_predicted_reduction << '\n';

      std::string scenario;
      bool accept = false;
      if (this->filter->acceptable(trial_progress.infeasibility, trial_merit)) {
         if (this->filter->acceptable_wrt_current_iterate(current_progress.infeasibility, current_merit, trial_progress.infeasibility, trial_merit)) {
            // switching condition: check whether the unconstrained predicted reduction is sufficiently positive
            if (this->switching_condition(merit_predicted_reduction, current_progress.infeasibility)) {
               // unconstrained Armijo sufficient decrease condition: predicted reduction should be positive (f-type)
               const double merit_actual_reduction = this->compute_actual_objective_reduction(current_merit, current_progress.infeasibility, trial_merit);
               DEBUG << "Unconstrained actual reduction = " << merit_actual_reduction << '\n';
               if (this->armijo_sufficient_decrease(merit_predicted_reduction, merit_actual_reduction)) {
                  DEBUG << "Trial iterate (f-type) was accepted by satisfying the Armijo condition\n";
                  accept = true;
               }
               else { // switching condition holds, but not Armijo condition
                  DEBUG << "Trial iterate (f-type) was rejected by violating the Armijo condition\n";
               }
               scenario = "f-type";
            }
               // switching condition violated: predicted reduction is not promising (h-type)
            else {
               DEBUG << "Trial iterate (h-type) was accepted by violating the switching condition\n";
               accept = true;
               this->filter->add(current_progress.infeasibility, current_merit);
               DEBUG << "Current iterate was added to the filter\n";
               scenario = "h-type";
            }
         }
         else {
            DEBUG << "Trial iterate not acceptable with respect to current point\n";
            scenario = "current";
         }
      }
      else {
         DEBUG << "Trial iterate not filter acceptable\n";
         scenario = "filter";
      }
      statistics.set("status", std::string(accept ? "✔" : "✘") + " (" + scenario + ")");
      return accept;
   }

   bool FletcherFilterMethod::is_infeasibility_sufficiently_reduced(const ProgressMeasures& /*reference_progress*/, const ProgressMeasures& trial_progress) const {
      // if the trial infeasibility improves upon the best known infeasibility
      return this->filter->infeasibility_sufficient_reduction(this->filter->get_smallest_infeasibility(), trial_progress.infeasibility);
   }
} // namespace
