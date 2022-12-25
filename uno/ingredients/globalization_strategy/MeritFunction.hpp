// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MERITFUNCTION_H
#define UNO_MERITFUNCTION_H

#include "GlobalizationStrategy.hpp"
#include "tools/Options.hpp"

class MeritFunction : public GlobalizationStrategy {
public:
   explicit MeritFunction(const Options& options);

   void initialize(const Iterate& first_iterate) override;
   [[nodiscard]] bool is_iterate_acceptable(const ProgressMeasures& current_progress, const ProgressMeasures& trial_progress,
         const PredictedReduction& predicted_reduction, double objective_multiplier) override;
   [[nodiscard]] bool is_feasibility_iterate_acceptable(double trial_infeasibility_measure) const override;
   void reset() override;
   void register_current_progress(const ProgressMeasures& current_progress) override;

protected:
   double smallest_known_infeasibility{INF<double>};
};

#endif // UNO_MERITFUNCTION_H
