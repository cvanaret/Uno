// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "WaechterFilterMethod.hpp"
#include "filters/Filter.hpp"
#include "ingredients/globalization_strategies/ProgressMeasures.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Logger.hpp"
#include "options/Options.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   WaechterFilterMethod::WaechterFilterMethod(const Options& options):
         FilterMethod(options),
         sufficient_infeasibility_decrease_factor(options.get_double("filter_sufficient_infeasibility_decrease_factor")) {
   }

   WaechterFilterMethod::~WaechterFilterMethod() { }

   void WaechterFilterMethod::initialize(Statistics& statistics, const Iterate& initial_iterate, const Options& options) {
      this->initial_infeasibility = initial_iterate.progress.infeasibility;
      FilterMethod::initialize(statistics, initial_iterate, options);
   }

   bool WaechterFilterMethod::is_regular_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
         const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction) {
      // in filter methods, we construct an unconstrained measure by ignoring infeasibility and scaling the objective measure by 1
      const double current_merit = FilterMethod::unconstrained_merit_function(current_progress);
      const double trial_merit = FilterMethod::unconstrained_merit_function(trial_progress);
      const double merit_predicted_reduction = FilterMethod::unconstrained_merit_function(predicted_reduction);
      DEBUG << "Current (infeasibility, objective + auxiliary) = (" << current_progress.infeasibility << ", " << current_merit << ")\n";
      DEBUG << "Trial   (infeasibility, objective + auxiliary) = (" << trial_progress.infeasibility << ", " << trial_merit << ")\n";
      DEBUG << "Current filter:\n" << *this->filter;
      DEBUG << "Unconstrained predicted reduction = " << merit_predicted_reduction << '\n';

      std::string scenario;
      bool accept = false;
      if (this->filter->acceptable(trial_progress.infeasibility, trial_merit)) {
         const double merit_actual_reduction = this->compute_actual_objective_reduction(current_merit, current_progress.infeasibility, trial_merit);
         DEBUG << "Unconstrained actual reduction = " << merit_actual_reduction << '\n';

         // TODO put this coefficient in the option file
         const bool small_infeasibility = current_progress.infeasibility <= 1e-4 * std::max(1., this->initial_infeasibility);
         const bool switching = (0. < merit_predicted_reduction) && this->switching_condition(merit_predicted_reduction, current_progress.infeasibility);
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
            scenario = "f-type";
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
         scenario = "filter";
      }
      statistics.set("status", std::string(accept ? "✔" : "✘") + " (" + scenario + ")");
      return accept;
   }

   bool WaechterFilterMethod::is_infeasibility_sufficiently_reduced(const ProgressMeasures& reference_progress, const ProgressMeasures& trial_progress) const {
      return trial_progress.infeasibility <= this->sufficient_infeasibility_decrease_factor * reference_progress.infeasibility &&
         this->filter->acceptable(trial_progress.infeasibility, FilterMethod::unconstrained_merit_function(trial_progress));
   }
} // namespace
