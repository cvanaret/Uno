// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "ingredients/globalization_strategies/ProgressMeasures.hpp"
#include "ingredients/globalization_strategies/switching_methods/filter_methods/FletcherFilterMethod.hpp"
#include "optimization/Iterate.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"
#include "tools/Statistics.hpp"

using namespace uno;

namespace {
   ProgressMeasures progress(double infeasibility, double objective) {
      return {infeasibility, [=](double) { return objective; }, 0., 1.};
   }

   Iterate trial_iterate(double infeasibility, double objective) {
      Iterate iterate(1, 1);
      iterate.progress = progress(infeasibility, objective);
      iterate.scaled_residuals.primal_infeasibility = infeasibility;
      return iterate;
   }
} // namespace

// Bug: the Fletcher filter left the restoration phase only if h < beta tau, where tau is the smallest infeasibility of the
// optimality filter. This test, stronger than FLT's requirement (a point acceptable to the filter), avoids switching
// phases too often. However, the filter may contain an entry whose infeasibility is at roundoff level (e.g. an h-type
// iterate at a point that is feasible up to roundoff): restoration then needs h ~ 1e-15 and may never hand back.
// Once the restoration iterate is essentially feasible (h <= filter_switch_back_infeasibility_threshold), filter
// acceptability (FLT's condition) suffices.
TEST(FletcherRestorationExit, EssentiallyFeasibleAcceptableIterateLeavesRestoration) {
   Options options;
   DefaultOptions::load(options);
   Presets::set(options, "filtersqp");
   FletcherFilterMethod filter_method(options);
   Statistics statistics{};
   statistics.add_column("Status", Statistics::string_width, 0);
   Iterate initial_iterate(1, 1);
   initial_iterate.progress = progress(1., 1.);
   filter_method.initialize(statistics, initial_iterate);
   // entries of the optimality filter: one at roundoff level, one entered before restoration
   filter_method.avoid_cycling_back_to(progress(1e-15, 1.));
   filter_method.avoid_cycling_back_to(progress(1e-3, 0.5));
   const double reference_infeasibility = 1e-3;

   // essentially feasible and acceptable to the filter: leaves restoration
   EXPECT_TRUE(filter_method.is_infeasibility_sufficiently_reduced(trial_iterate(1e-10, 0.9), reference_infeasibility));
   // essentially feasible but dominated by the roundoff entry: stays in restoration
   EXPECT_FALSE(filter_method.is_infeasibility_sufficiently_reduced(trial_iterate(1e-10, 2.), reference_infeasibility));
   // not essentially feasible: the strict test h < beta tau applies
   EXPECT_FALSE(filter_method.is_infeasibility_sufficiently_reduced(trial_iterate(1e-5, 0.1), reference_infeasibility));
}

// the strict test is unchanged when the smallest infeasibility of the filter is not at roundoff level
TEST(FletcherRestorationExit, StrictTestWithoutRoundoffEntries) {
   Options options;
   DefaultOptions::load(options);
   Presets::set(options, "filtersqp");
   FletcherFilterMethod filter_method(options);
   Statistics statistics{};
   statistics.add_column("Status", Statistics::string_width, 0);
   Iterate initial_iterate(1, 1);
   initial_iterate.progress = progress(1., 1.);
   filter_method.initialize(statistics, initial_iterate);
   filter_method.avoid_cycling_back_to(progress(1e-3, 0.5));
   EXPECT_TRUE(filter_method.is_infeasibility_sufficiently_reduced(trial_iterate(5e-4, 10.), 1e-3));
   EXPECT_FALSE(filter_method.is_infeasibility_sufficiently_reduced(trial_iterate(1e-3, 0.1), 1e-3));
}
