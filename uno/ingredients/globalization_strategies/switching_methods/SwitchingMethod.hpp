// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SWITCHINGMETHOD_H
#define UNO_SWITCHINGMETHOD_H

#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"

namespace uno {
   class SwitchingMethod : public GlobalizationStrategy {
   public:
      explicit SwitchingMethod(const Options& options);
      ~SwitchingMethod() override = default;

   protected:
      const double delta;
      const double switching_merit_exponent;
      const double switching_infeasibility_exponent;

      [[nodiscard]] static double unconstrained_merit_function(const ProgressMeasures& progress);
      [[nodiscard]] bool switching_condition(double predicted_reduction, double step_length, double current_infeasibility) const;
      [[nodiscard]] double compute_actual_merit_reduction(double current_merit, double trial_merit) const;
   };
} // namespace

#endif // UNO_SWITCHINGMETHOD_H