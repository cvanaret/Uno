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
   bool is_iterate_acceptable(const ProgressMeasures& current_progress, const ProgressMeasures& trial_progress, const ProgressMeasures& predicted_reduction) override;
   void reset() override;
   void register_current_progress(const ProgressMeasures& current_progress) override;
};

#endif // UNO_MERITFUNCTION_H
