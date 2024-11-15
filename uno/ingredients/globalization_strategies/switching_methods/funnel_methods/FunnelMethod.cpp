// Copyright (c) 2024 David Kiessling, Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "FunnelMethod.hpp"
#include "../../ProgressMeasures.hpp"
#include "optimization/Iterate.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   FunnelMethod::FunnelMethod(const Options& options) :
         SwitchingMethod(options),
         funnel(options),
         parameters({
               options.get_double("funnel_ubd"),
               options.get_double("funnel_fact"),
               options.get_double("funnel_beta"),
               options.get_double("funnel_gamma"),
         }),
         require_acceptance_wrt_current_iterate(options.get_bool("funnel_require_acceptance_wrt_current_iterate")) {
   }

   void FunnelMethod::initialize(Statistics& statistics, const Iterate& initial_iterate, const Options& options) {
      const double upper_bound = std::max(this->parameters.initial_upper_bound, this->parameters.infeasibility_factor * initial_iterate.progress.infeasibility);
      this->funnel.set_infeasibility_upper_bound(upper_bound);
      DEBUG << "Current funnel width: " << this->funnel.current_width() << '\n';

      statistics.add_column("funnel width", Statistics::double_width, options.get_int("statistics_funnel_width_column_order"));
      statistics.set("funnel width", this->funnel.current_width());
   }

   bool FunnelMethod::is_regular_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
         const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction) {
      bool accept = false;
      std::string scenario;
      if (this->funnel.acceptable(trial_progress.infeasibility)) {
         // in filter and funnel methods, we construct an unconstrained measure by ignoring infeasibility and scaling the objective measure by 1
         const double current_merit = SwitchingMethod::unconstrained_merit_function(current_progress);
         const double trial_merit = SwitchingMethod::unconstrained_merit_function(trial_progress);
         const double merit_predicted_reduction = SwitchingMethod::unconstrained_merit_function(predicted_reduction);
         DEBUG << "Current: (infeasibility, objective + auxiliary) = (" << current_progress.infeasibility << ", " << current_merit << ")\n";
         DEBUG << "Trial:   (infeasibility, objective + auxiliary) = (" << trial_progress.infeasibility << ", " << trial_merit << ")\n";
         DEBUG << "Unconstrained predicted reduction = " << merit_predicted_reduction << '\n';

         // IF require_acceptance_wrt_current_iterate == false, then condition always fulfilled, we never check
         // If true, then first part is false, and we always check wrt current iterate
         if (not this->require_acceptance_wrt_current_iterate ||
             this->acceptable_wrt_current_iterate(current_progress.infeasibility, current_merit, trial_progress.infeasibility, trial_merit)) {
            // f-type step
            if (this->switching_condition(merit_predicted_reduction, current_progress.infeasibility)) {
               DEBUG << "Trial iterate satisfies switching condition\n";
               // unconstrained Armijo sufficient decrease condition (predicted reduction should be positive)
               const double objective_actual_reduction = this->compute_actual_objective_reduction(current_merit, trial_merit);
               DEBUG << "Unconstrained actual reduction = " << objective_actual_reduction << '\n';
               if (this->armijo_sufficient_decrease(merit_predicted_reduction, objective_actual_reduction)) {
                  DEBUG << "Trial iterate (f-type) was ACCEPTED by satisfying Armijo condition\n";
                  accept = true;
               }
               else { // switching condition holds, but not Armijo condition
                  DEBUG << "Trial iterate (f-type) was REJECTED by violating the Armijo condition\n";
               }
               scenario = "f-type";
            }
            // h-type step
            else if (this->funnel.sufficient_decrease_condition(trial_progress.infeasibility)) {
               DEBUG << "\t\tTrial iterate  (h-type) ACCEPTED by violating the switching condition ...\n";
               accept = true;

               DEBUG << "\t\tEntering funnel reduction mechanism\n";
               this->funnel.update(current_progress.infeasibility, trial_progress.infeasibility);
               statistics.set("funnel width", this->funnel.current_width());
               scenario = "h-type";
            }
            else {
               DEBUG << "\t\tTrial iterate REJECTED by violating switching and funnel sufficient decrease condition\n";
               scenario = "funnel";
            }
         }
         else {
            DEBUG << "Trial iterate not acceptable with respect to current point\n";
            scenario = "current";
         }
      }
      else {
         DEBUG << "\t\tTrial iterate REJECTED. Not in funnel\n";
         scenario = "funnel";
      }

      statistics.set("status", std::string(accept ? "✔" : "✘") + " (" + scenario + ")");
      if (accept) {
         this->funnel.print();
      }
      return accept;
   }

   bool FunnelMethod::is_infeasibility_sufficiently_reduced(const ProgressMeasures& reference_progress, const ProgressMeasures& trial_progress) const {
      return this->funnel.acceptable(trial_progress.infeasibility) &&
             trial_progress.infeasibility <= this->parameters.beta * reference_progress.infeasibility;
   }

   void FunnelMethod::reset() {
      // do nothing
   }

   void FunnelMethod::notify_switch_to_feasibility(const ProgressMeasures& /*current_progress_measures*/) {
   }

   void FunnelMethod::notify_switch_to_optimality(const ProgressMeasures& current_progress_measures) {
      DEBUG << "Funnel is reduced after restoration phase\n";
      this->funnel.update_restoration(current_progress_measures.infeasibility);
   }

   // check acceptability wrt current point
   bool FunnelMethod::acceptable_wrt_current_iterate(double current_infeasibility, double current_objective, double trial_infeasibility,
         double trial_objective) const {
      return this->infeasibility_sufficient_reduction(current_infeasibility, trial_infeasibility) ||
         this->objective_sufficient_reduction(current_objective, trial_objective, trial_infeasibility);
   }

   bool FunnelMethod::infeasibility_sufficient_reduction(double current_infeasibility, double trial_infeasibility) const {
      return (trial_infeasibility < this->parameters.beta * current_infeasibility);
   }

   bool FunnelMethod::objective_sufficient_reduction(double current_objective, double trial_objective, double trial_infeasibility) const {
      return (trial_objective <= current_objective - this->parameters.gamma * trial_infeasibility);
   }

   double FunnelMethod::compute_actual_objective_reduction(double current_objective_measure, double trial_objective_measure) {
      double actual_reduction = current_objective_measure - trial_objective_measure;
      if (this->protect_actual_reduction_against_roundoff) {
         static double machine_epsilon = std::numeric_limits<double>::epsilon();
         actual_reduction += 10. * machine_epsilon * std::abs(current_objective_measure);
      }
      return actual_reduction;
   }

   void FunnelMethod::set_statistics(Statistics& statistics) const {
      statistics.set("funnel width", this->funnel.current_width());
   }
} // namespace
