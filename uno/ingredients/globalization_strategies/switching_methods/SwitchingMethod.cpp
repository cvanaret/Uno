// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include "SwitchingMethod.hpp"
#include "../ProgressMeasures.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Logger.hpp"
#include "options/Options.hpp"

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
         const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction, double /*objective_multiplier*/) {
      this->set_statistics(statistics);
      return this->is_regular_iterate_acceptable(statistics, current_progress, trial_progress, predicted_reduction);
   }
} // namespace