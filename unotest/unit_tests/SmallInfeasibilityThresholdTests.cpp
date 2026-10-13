// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// theta_min (IPOPT): theta_min_fact * max(1, theta(x_0)), shared by the switching rule and the minimum step length
// of the line search (Wächter & Biegler 2006, eq. (23)). x_0 is the point at which the filter is initialized.

#include <gtest/gtest.h>
#include "ingredients/globalization_strategies/switching_methods/filter_methods/WaechterFilterMethod.hpp"
#include "optimization/Iterate.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"
#include "tools/Statistics.hpp"

using namespace uno;

namespace {
   Options ipopt_options() {
      Options options;
      DefaultOptions::load(options);
      Presets::set(options, "ipopt");
      return options; // barrier_small_infeasibility_factor = 1e-4
   }

   double threshold_after_initialization(WaechterFilterMethod& filter_method, double initial_infeasibility) {
      Statistics statistics;
      Iterate iterate(1, 1);
      iterate.progress.infeasibility = initial_infeasibility;
      filter_method.initialize(statistics, iterate);
      return filter_method.small_infeasibility_threshold();
   }
} // namespace

TEST(SmallInfeasibilityThreshold, RelativeToInitialInfeasibility) {
   const Options options = ipopt_options();
   WaechterFilterMethod filter_method(options);
   EXPECT_DOUBLE_EQ(threshold_after_initialization(filter_method, 1e3), 1e-1);
}

TEST(SmallInfeasibilityThreshold, AbsoluteBelowOne) {
   const Options options = ipopt_options();
   WaechterFilterMethod filter_method(options);
   EXPECT_DOUBLE_EQ(threshold_after_initialization(filter_method, 0.5), 1e-4);
}

// the restoration filter is initialized at the first restoration iterate: its theta_min follows
TEST(SmallInfeasibilityThreshold, RecomputedUponReinitialization) {
   const Options options = ipopt_options();
   WaechterFilterMethod filter_method(options);
   threshold_after_initialization(filter_method, 1e3);
   EXPECT_DOUBLE_EQ(threshold_after_initialization(filter_method, 50.), 5e-3);
}

// IPOPT does not recompute theta_min when the filter is reset (FilterLSAcceptor::Reset only clears the filter)
TEST(SmallInfeasibilityThreshold, UnchangedByFilterReset) {
   const Options options = ipopt_options();
   WaechterFilterMethod filter_method(options);
   const double threshold = threshold_after_initialization(filter_method, 1e3);
   filter_method.reset();
   EXPECT_DOUBLE_EQ(filter_method.small_infeasibility_threshold(), threshold);
}
