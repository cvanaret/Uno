// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Bug: the Waechter filter method accepted h-type trial iterates that reduce the infeasibility while increasing the
// merit by an arbitrary amount (CUTEst LUKVLI10: from 1e4 to 1e61 in one step, after which the iterates diverged).
// IPOPT rejects them (option obj_max_inc): the merit may not increase by more than filter_merit_max_increase orders of
// magnitude (relative to its current magnitude)

#include <gtest/gtest.h>
#include "ingredients/globalization_strategies/ProgressMeasures.hpp"
#include "ingredients/globalization_strategies/switching_methods/filter_methods/WaechterFilterMethod.hpp"
#include "optimization/Iterate.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "tools/Statistics.hpp"

using namespace uno;

namespace {
   ProgressMeasures progress(double infeasibility, double objective) {
      return {infeasibility, [=](double) { return objective; }, 0., 1.};
   }

   Decision test_h_type_step(double current_objective, double trial_objective) {
      Options options;
      DefaultOptions::load(options);
      Statistics statistics{};
      statistics.add_column("Status", Statistics::string_width, 0);
      WaechterFilterMethod method(options);
      Iterate initial_iterate(1, 1);
      initial_iterate.progress = progress(10., current_objective);
      method.initialize(statistics, initial_iterate);
      // the infeasibility is halved; zero predicted reduction: the switching condition does not hold (h-type)
      return method.is_iterate_acceptable(statistics, progress(10., current_objective), progress(5., trial_objective),
         progress(0., 0.), 1.);
   }
} // namespace

TEST(FilterMeritIncrease, ModerateIncreaseIsAccepted) {
   EXPECT_TRUE(test_h_type_step(1e4, 1e5).accepted);
   EXPECT_TRUE(test_h_type_step(0., 1e5).accepted);
}

TEST(FilterMeritIncrease, AstronomicalIncreaseIsRejected) {
   // 1e4 -> 1e61 (LUKVLI10): log10(increase) = 61 > 5 + log10(1e4)
   EXPECT_FALSE(test_h_type_step(1e4, 1e61).accepted);
   // below magnitude 10, the reference magnitude is 1: increases beyond 1e6 are rejected
   EXPECT_FALSE(test_h_type_step(0., 1e7).accepted);
}
