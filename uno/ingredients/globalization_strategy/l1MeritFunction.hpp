// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MERITFUNCTION_H
#define UNO_MERITFUNCTION_H

#include "GlobalizationStrategy.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   class l1MeritFunction : public GlobalizationStrategy {
   public:
      explicit l1MeritFunction(const Options& options);

      void initialize(Statistics& statistics, const Iterate& initial_iterate, const Options& options) override;
      [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
            const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction, double objective_multiplier) override;
      [[nodiscard]] bool is_infeasibility_sufficiently_reduced(const ProgressMeasures& current_progress, const ProgressMeasures& trial_progress) const override;
      void reset() override;
      void notify_switch_to_feasibility(const ProgressMeasures& current_progress) override;
      void notify_switch_to_optimality(const ProgressMeasures& current_progress) override;

   protected:
      double smallest_known_infeasibility{INF<double>};

      [[nodiscard]] static double constrained_merit_function(const ProgressMeasures& progress, double objective_multiplier);
      [[nodiscard]] double compute_merit_actual_reduction(double current_merit_value, double trial_merit_value) const;
   };
} // namespace

#endif // UNO_MERITFUNCTION_H
