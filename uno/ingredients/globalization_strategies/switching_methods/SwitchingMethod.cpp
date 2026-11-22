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

   // an h-type iterate makes progress towards feasibility, which is impossible at a feasible iterate. There, the
   // switching condition reduces to pred > 0, and the trial iterate must be f-type (Armijo). If pred <= 0 (e.g. a KKT
   // point of a nonconvex QP with a large trust region, or a non-descent line-search direction), the trial iterate is
   // rejected: the trust region shrinks until the QP step predicts a decrease of order radius * stationarity measure
   bool SwitchingMethod::allows_h_type_iterate(double current_infeasibility) {
      return (0. < current_infeasibility);
   }
} // namespace
