// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "WaechterFilterMethod.hpp"
#include "filters/Filter.hpp"
#include "ingredients/globalization_strategies/ProgressMeasures.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Logger.hpp"
#include "options/Options.hpp"
#include "tools/Statistics.hpp"
#include "tools/Symbols.hpp"

namespace uno {
   WaechterFilterMethod::WaechterFilterMethod(const Options& options):
         FilterMethod(options),
         sufficient_infeasibility_decrease_factor(options.get_double("filter_sufficient_infeasibility_decrease_factor")),
         small_infeasibility_factor(options.get_double("barrier_small_infeasibility_factor")),
         filter_reset_iteration_threshold(options.get_unsigned_int("filter_reset_iteration_threshold")),
         max_number_filter_resets(options.get_unsigned_int("max_number_filter_resets")) {
   }

   void WaechterFilterMethod::initialize(Statistics& statistics, const Iterate& initial_iterate) {
      this->initial_infeasibility = initial_iterate.progress.infeasibility;
      FilterMethod::initialize(statistics, initial_iterate);
   }

   bool WaechterFilterMethod::is_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
         const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction, double /*objective_multiplier*/) {
      // in filter methods, we construct an unconstrained measure by ignoring infeasibility and scaling the objective measure by 1
      const double current_merit = FilterMethod::unconstrained_merit_function(current_progress);
      const double trial_merit = FilterMethod::unconstrained_merit_function(trial_progress);
      const double merit_predicted_reduction = FilterMethod::unconstrained_merit_function(predicted_reduction);
      DEBUG << "Current (infeasibility, objective + auxiliary) = (" << current_progress.infeasibility << ", " << current_merit << ")\n";
      DEBUG << "Trial   (infeasibility, objective + auxiliary) = (" << trial_progress.infeasibility << ", " << trial_merit << ")\n";
      DEBUG << "Current filter:\n" << *this->filter;
      DEBUG << "Unconstrained predicted reduction = " << merit_predicted_reduction << '\n';

      if (!this->filter->acceptable_wrt_infeasibility_upper_bound(trial_progress.infeasibility)) {
         DEBUG << "Trial iterate not acceptable wrt infeasibility upper bound\n";
         statistics.set("Status", std::string(symbols::fail) + " (upper bound)");
         return false;
      }
      // now acceptable wrt infeasibility upper bound

      const double merit_actual_reduction = this->compute_actual_objective_reduction(current_merit, current_progress.infeasibility,
         trial_merit);
      DEBUG << "Unconstrained actual reduction = " << merit_actual_reduction << '\n';

      const bool small_infeasibility = current_progress.infeasibility <= this->small_infeasibility_factor *
         std::max(1., this->initial_infeasibility);
      const bool switching = (0. < merit_predicted_reduction) && this->switching_condition(merit_predicted_reduction,
         current_progress.infeasibility);
      const bool sufficient_decrease = this->armijo_sufficient_decrease(merit_predicted_reduction, merit_actual_reduction);

      // switching condition: the unconstrained predicted reduction is sufficiently positive
      if (switching && small_infeasibility) {
         DEBUG << "Switching condition satisfied\n";
         // unconstrained Armijo sufficient decrease condition: predicted reduction should be positive (f-type)
         if (sufficient_decrease) {
            DEBUG << "Trial iterate (f-type) was accepted by satisfying Armijo condition\n";
            statistics.set("Status", std::string(symbols::check) + " (f-type)");
         }
         else {
            DEBUG << "Armijo condition not satisfied\n";
            statistics.set("Status", std::string(symbols::fail) + " (Armijo)");
            this->last_rejection_due_to_filter = false;
            return false;
         }
      }
      else {
         DEBUG << "Switching condition violated\n";
         if (this->filter->acceptable_wrt_current_iterate(current_progress.infeasibility, current_merit,
               trial_progress.infeasibility, trial_merit)) {
            DEBUG << "Trial iterate (h-type) acceptable with respect to current point\n";
            statistics.set("Status", std::string(symbols::check) + " (h-type)");
         }
         else {
            DEBUG << "Trial iterate (h-type) not acceptable with respect to current point\n";
            statistics.set("Status", std::string(symbols::fail) + " (current)");
            this->last_rejection_due_to_filter = false;
            return false;
         }
      }

      // check acceptability wrt filter
      if (!this->filter->filter_acceptable(trial_progress.infeasibility, trial_merit)) {
         DEBUG << "Trial iterate not filter acceptable\n";
         statistics.set("Status", std::string(symbols::fail) + " (filter)");
         this->last_rejection_due_to_filter = true;
         return false;
      }
      // now filter acceptable

      // filter reset heuristic
      if (this->number_filter_resets < this->max_number_filter_resets) {
         if (this->last_rejection_due_to_filter) {
            ++this->number_successive_filter_rejections;
            if (this->number_successive_filter_rejections >= this->filter_reset_iteration_threshold) {
               DEBUG << "Resetting the filter\n";
               this->filter->reset();
               ++this->number_filter_resets;
            }
         }
         else {
            this->number_successive_filter_rejections = 0;
         }
      }
      else {
         DEBUG << "The filter cannot be reset (max number of resets exceeded)\n";
      }
      this->last_rejection_due_to_filter = false;

      // possibly augment the filter
      if (!switching || !sufficient_decrease) {
         DEBUG << "Adding current iterate to the filter\n";
         this->filter->add(current_progress.infeasibility, current_merit);
      }
      return true;
   }

   bool WaechterFilterMethod::is_infeasibility_sufficiently_reduced(const Iterate& trial_iterate, double reference_infeasibility) const {
      // use the infeasibility in the residual norm here (inf norm for IPOPT implementation)
      DEBUG << "Testing whether " << trial_iterate.primal_infeasibility << " <= " << this->sufficient_infeasibility_decrease_factor
         << "*" << reference_infeasibility << " = " << this->sufficient_infeasibility_decrease_factor * reference_infeasibility << '\n';
      return trial_iterate.primal_infeasibility <= this->sufficient_infeasibility_decrease_factor * reference_infeasibility &&
         this->filter->filter_acceptable(trial_iterate.progress.infeasibility, unconstrained_merit_function(trial_iterate.progress));
   }

   std::string WaechterFilterMethod::get_name() const {
      return "Waechter-filter";
   }
} // namespace
