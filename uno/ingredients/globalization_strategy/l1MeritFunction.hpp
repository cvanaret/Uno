// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MERITFUNCTION_H
#define UNO_MERITFUNCTION_H

#include "GlobalizationStrategy.hpp"

class l1MeritFunction : public GlobalizationStrategy {
public:
   explicit l1MeritFunction(Statistics& statistics, const Options& options);

   void initialize(const Iterate& initial_iterate) override;
   [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, const Iterate& trial_iterate, const ProgressMeasures& current_progress,
         const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction, double objective_multiplier) override;
   [[nodiscard]] bool is_infeasibility_acceptable(double infeasibility_measure) const override;
   void reset() override;
   void register_current_progress(const ProgressMeasures& current_progress) override;

protected:
   double smallest_known_infeasibility{INF<double>};
};

#endif // UNO_MERITFUNCTION_H
