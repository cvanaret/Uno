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

      bool is_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress, const ProgressMeasures& trial_progress,
            const ProgressMeasures& predicted_reduction, double objective_multiplier) override;

   protected:
      const double delta;
      const double switching_infeasibility_exponent;

      [[nodiscard]] virtual bool is_regular_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
            const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction) = 0;
      [[nodiscard]] static double unconstrained_merit_function(const ProgressMeasures& progress);
      [[nodiscard]] bool switching_condition(double predicted_reduction, double current_infeasibility) const;
      [[nodiscard]] bool is_feasibility_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
            const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction);
      virtual void set_statistics(Statistics& statistics) const = 0;
   };
} // namespace

#endif // UNO_SWITCHINGMETHOD_H