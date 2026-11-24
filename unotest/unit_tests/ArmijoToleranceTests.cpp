// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Wächter & Biegler's f-type condition (11)/(15) is phi(trial) <= phi(current) + eta * m(alpha), with eta = 1e-8 in the
// IPOPT preset. The absolute armijo_tolerance (1e-9) is subtracted from the predicted reduction: below it, the Armijo
// condition degenerates into simple decrease, and the filter is not augmented. Whether that happens depends on the
// scaling of the problem, not on the progress of the method. IPOPT has no such tolerance (only a relative roundoff
// guard on the actual reduction, protect_actual_reduction_against_roundoff in Uno)

#include <gtest/gtest.h>
#include "ingredients/globalization_strategies/ProgressMeasures.hpp"
#include "ingredients/globalization_strategies/switching_methods/filter_methods/WaechterFilterMethod.hpp"
#include "optimization/Iterate.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"
#include "tools/Statistics.hpp"

using namespace uno;

namespace {
   ProgressMeasures progress(double infeasibility, double objective, double step_length = 1.) {
      return {infeasibility, [=](double) { return objective; }, 0., step_length};
   }

   Options ipopt_options() {
      Options options;
      DefaultOptions::load(options);
      Presets::set(options, "ipopt");
      return options;
   }

   // f-type step at a feasible iterate, for a problem whose objective is of order 1e-6 (e.g. after gradient-based
   // scaling) and whose predicted reduction is 5e-10
   bool f_type_step_is_accepted(double actual_reduction) {
      const Options options = ipopt_options();
      WaechterFilterMethod strategy(options);
      Statistics statistics{};
      statistics.add_column("Status", Statistics::string_width, 0);
      constexpr double objective = 1e-6;
      Iterate initial_iterate(1, 1);
      initial_iterate.progress = progress(0., objective);
      strategy.initialize(statistics, initial_iterate);
      constexpr double predicted_reduction = 5e-10;
      return strategy.is_iterate_acceptable(statistics, progress(0., objective), progress(0., objective - actual_reduction),
         progress(0., predicted_reduction), 1.).accepted;
   }
} // namespace

// the trial iterate does not decrease the objective: the Armijo condition (eta * pred = 5e-18 > 0) does not hold
TEST(ArmijoTolerance, IpoptPresetRejectsStepWithoutDecrease) {
   EXPECT_FALSE(f_type_step_is_accepted(0.));
}

// a decrease of 1e-20 is far below eta * pred = 5e-18 (and above the roundoff guard 10 eps |phi| ~ 2e-21)
TEST(ArmijoTolerance, IpoptPresetRejectsInsufficientDecrease) {
   EXPECT_FALSE(f_type_step_is_accepted(1e-20));
}
