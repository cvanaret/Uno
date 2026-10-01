// Copyright (c) 2024-2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include "SwitchingMethod.hpp"
#include "../ProgressMeasures.hpp"
#include "options/Options.hpp"

namespace uno {
   SwitchingMethod::SwitchingMethod(const Options& options): GlobalizationStrategy(options),
      delta(options.get_double("switching_delta")),
      switching_merit_exponent(options.get_double("switching_merit_exponent")),
      switching_infeasibility_exponent(options.get_double("switching_infeasibility_exponent")) { }

   double SwitchingMethod::unconstrained_merit_function(const ProgressMeasures& progress) {
      return progress.objective(1.) + progress.auxiliary;
   }

   // IPOPT, eq. (19): α (-∇φᵀd)^{sφ} > δ θ^{sθ}
   // in Uno, m(α) = -α ∇φᵀd, therefore the LHS is given by α (m(α)/α)^{sφ} > δ θ^{sθ}
   bool SwitchingMethod::switching_condition(double predicted_reduction, double step_length, double current_infeasibility) const {
      if (predicted_reduction <= 0.) {
         return false;
      }
      const double unit_step_reduction = predicted_reduction / step_length;
      return step_length * std::pow(unit_step_reduction, this->switching_merit_exponent) >
         this->delta * std::pow(current_infeasibility, this->switching_infeasibility_exponent);
   }

   double SwitchingMethod::compute_actual_merit_reduction(double current_merit, double trial_merit) const {
      double actual_reduction = current_merit - trial_merit;
      if (this->protect_actual_reduction_against_roundoff) {
         static double machine_epsilon = std::numeric_limits<double>::epsilon();
         actual_reduction += this->protected_actual_reduction_macheps_coefficient * machine_epsilon * std::abs(current_merit);
      }
      return actual_reduction;
   }
} // namespace