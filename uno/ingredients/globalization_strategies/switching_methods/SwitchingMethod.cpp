// Copyright (c) 2024-2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include "SwitchingMethod.hpp"
#include "../ProgressMeasures.hpp"
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
} // namespace