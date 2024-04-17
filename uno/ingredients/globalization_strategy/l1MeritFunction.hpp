// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MERITFUNCTION_H
#define UNO_MERITFUNCTION_H

#include "GlobalizationStrategy.hpp"

class l1MeritFunction : public GlobalizationStrategy {
public:
   explicit l1MeritFunction(const Options& options);

   void initialize(Statistics& statistics, const Iterate& initial_iterate, const Options& options) override;
   [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, const ProgressMeasures& current_progress,
         const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction, double objective_multiplier) override;
   [[nodiscard]] bool is_infeasibility_acceptable(const ProgressMeasures& current_progress, const ProgressMeasures& trial_progress) const override;
   void reset() override;
   void register_current_progress(const ProgressMeasures& current_progress) override;
   [[nodiscard]] double get_infeasibility_upper_bound() const override;
   void set_infeasibility_upper_bound(double new_upper_bound) override;

protected:
   double smallest_known_infeasibility{INF<double>};
};

#endif // UNO_MERITFUNCTION_H
