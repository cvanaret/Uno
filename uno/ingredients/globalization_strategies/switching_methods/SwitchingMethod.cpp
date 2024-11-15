// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include "SwitchingMethod.hpp"
#include "../ProgressMeasures.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Logger.hpp"
#include "options/Options.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   SwitchingMethod::SwitchingMethod(const Options& options): GlobalizationStrategy(options),
      delta(options.get_double("switching_delta")),
      switching_infeasibility_exponent(options.get_double("switching_infeasibility_exponent")) { }

   double SwitchingMethod::unconstrained_merit_function(const ProgressMeasures& progress) {
      return progress.objective(1.) + progress.auxiliary;
   }

   bool SwitchingMethod::switching_condition(double predicted_reduction, double current_infeasibility) const {
      return predicted_reduction > this->delta * std::pow(current_infeasibility, this->switching_infeasibility_exponent);
   }

   /* check acceptability of step
    * switching methods enforce an *unconstrained* sufficient decrease condition
    * precondition: feasible step
    * */
   bool SwitchingMethod::is_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
         const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction, double objective_multiplier) {
      this->set_statistics(statistics);
      const bool solving_feasibility_problem = (objective_multiplier == 0.);
      if (solving_feasibility_problem) {
         return this->is_feasibility_iterate_acceptable(statistics, current_progress, trial_progress, predicted_reduction);
      }
      else {
         return this->is_regular_iterate_acceptable(statistics, current_progress, trial_progress, predicted_reduction);
      }
   }

   // solving the feasibility problem = working on infeasibility only (no filter acceptability test)
   bool SwitchingMethod::is_feasibility_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
         const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction) {
      // drop the objective measure and focus on infeasibility and auxiliary terms (barrier, proximal, ...)
      const double current_merit = current_progress.infeasibility + current_progress.auxiliary;
      const double trial_merit = trial_progress.infeasibility + trial_progress.auxiliary;
      const double predicted_merit_reduction = predicted_reduction.infeasibility + predicted_reduction.auxiliary;
      const double actual_merit_reduction = current_merit - trial_merit;
      DEBUG << "Current merit = " << current_merit << '\n';
      DEBUG << "Trial merit = " << trial_merit << '\n';
      DEBUG << "Predicted merit reduction = " << predicted_merit_reduction << '\n';
      DEBUG << "Actual merit reduction = " << actual_merit_reduction << '\n';

      bool accept = false;
      if (this->armijo_sufficient_decrease(predicted_merit_reduction, actual_merit_reduction)) {
         DEBUG << "Trial iterate (h-type) was accepted by satisfying the Armijo condition\n";
         accept = true;
      }
      else {
         DEBUG << "Trial iterate (h-type) was rejected by violating the Armijo condition\n";
      }
      Iterate::number_eval_objective--;
      statistics.set("status", std::string(accept ? "✔" : "✘") + " (restoration)");
      return accept;
   }
} // namespace