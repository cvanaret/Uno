// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Bug: at a feasible iterate (h = 0), a feasible trial iterate (h = 0) with a nonpositive predicted reduction was accepted
// as an h-type iterate, whatever its objective. The switching condition fails (pred <= delta * 0 = 0), the trial
// iterate is "acceptable wrt the current iterate" because 0 <= beta * 0, and it is filter acceptable as the left-most
// entry (same non-strict test). The objective could therefore increase at every iteration (ascent, cycling).
// The FLT convergence theory excludes this case because the QP is solved globally: d = 0 is feasible at a feasible
// iterate, hence pred >= 0. Active-set QP solvers (BQPD) on nonconvex QPs only return KKT points, which can have
// pred <= 0 for large trust-region radii; so can a line-search direction without inertia correction.
// filterSQP is not affected: the infeasibility test is strict (hc < beta * cons), so a feasible trial iterate must
// decrease the objective.

#include <string>
#include <gtest/gtest.h>
#include "ingredients/globalization_strategies/GlobalizationStrategy.hpp"
#include "ingredients/globalization_strategies/ProgressMeasures.hpp"
#include "ingredients/globalization_strategies/switching_methods/filter_methods/Filter.hpp"
#include "ingredients/globalization_strategies/switching_methods/filter_methods/FletcherFilterMethod.hpp"
#include "ingredients/globalization_strategies/switching_methods/filter_methods/WaechterFilterMethod.hpp"
#include "ingredients/globalization_strategies/switching_methods/funnel_methods/FunnelMethod.hpp"
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

   Options preset_options(const std::string& preset) {
      Options options;
      DefaultOptions::load(options);
      Presets::set(options, preset);
      return options;
   }

   // drive the strategy along feasible iterates: each trial iterate is feasible, has a larger objective than the
   // current iterate and a nonpositive predicted reduction. Return the number of accepted (ascent) iterates
   template <typename Strategy>
   size_t number_accepted_ascent_steps(const std::string& preset, double predicted_reduction) {
      const Options options = preset_options(preset);
      Strategy strategy(options);
      Statistics statistics{};
      statistics.add_column("Status", Statistics::string_width, 0);
      statistics.add_column("Funnel", Statistics::double_width, 2);
      double objective = 1.;
      Iterate initial_iterate(1, 1);
      initial_iterate.progress = progress(0., objective);
      strategy.initialize(statistics, initial_iterate);

      size_t number_accepted = 0;
      for (size_t iteration = 0; iteration < 5; ++iteration) {
         const Decision decision = strategy.is_iterate_acceptable(statistics, progress(0., objective),
            progress(0., objective + 1.), progress(0., predicted_reduction), 1.);
         if (decision.accepted) {
            ++number_accepted;
            objective += 1.;
         }
      }
      return number_accepted;
   }

   // a feasible trial iterate that decreases the objective, with a positive predicted reduction (f-type)
   template <typename Strategy>
   bool descent_step_is_accepted(const std::string& preset) {
      const Options options = preset_options(preset);
      Strategy strategy(options);
      Statistics statistics{};
      statistics.add_column("Status", Statistics::string_width, 0);
      statistics.add_column("Funnel", Statistics::double_width, 2);
      Iterate initial_iterate(1, 1);
      initial_iterate.progress = progress(0., 1.);
      strategy.initialize(statistics, initial_iterate);
      return strategy.is_iterate_acceptable(statistics, progress(0., 1.), progress(0., 0.5), progress(0., 0.5), 1.).accepted;
   }
} // namespace

TEST(FeasiblePointAcceptance, FilterRequiresObjectiveDecreaseAtZeroInfeasibility) {
   Filter filter(preset_options("filtersqp"));
   filter.set_infeasibility_upper_bound(1e2);
   filter.add(0., 1.);
   // a feasible entry dominates every feasible point with a larger objective
   EXPECT_FALSE(filter.filter_acceptable(0., 2.));
   EXPECT_FALSE(filter.acceptable_wrt_current_iterate(0., 1., 0., 2.));
   EXPECT_TRUE(filter.filter_acceptable(0., 0.5));
   EXPECT_TRUE(filter.acceptable_wrt_current_iterate(0., 1., 0., 0.5));
   // positive infeasibilities are unaffected
   EXPECT_TRUE(filter.acceptable_wrt_current_iterate(1., 1., 0.5, 2.));
}

TEST(FeasiblePointAcceptance, FletcherFilterRejectsAscentAtFeasiblePoint) {
   EXPECT_EQ(number_accepted_ascent_steps<FletcherFilterMethod>("filtersqp", -0.1), 0);
   EXPECT_EQ(number_accepted_ascent_steps<FletcherFilterMethod>("filtersqp", 0.), 0);
   EXPECT_TRUE(descent_step_is_accepted<FletcherFilterMethod>("filtersqp"));
}

TEST(FeasiblePointAcceptance, WaechterFilterRejectsAscentAtFeasiblePoint) {
   EXPECT_EQ(number_accepted_ascent_steps<WaechterFilterMethod>("ipopt", -0.1), 0);
   EXPECT_EQ(number_accepted_ascent_steps<WaechterFilterMethod>("ipopt", 0.), 0);
   EXPECT_TRUE(descent_step_is_accepted<WaechterFilterMethod>("ipopt"));
}

TEST(FeasiblePointAcceptance, FunnelRejectsAscentAtFeasiblePoint) {
   EXPECT_EQ(number_accepted_ascent_steps<FunnelMethod>("funnelsqp", -0.1), 0);
   EXPECT_EQ(number_accepted_ascent_steps<FunnelMethod>("funnelsqp", 0.), 0);
   EXPECT_TRUE(descent_step_is_accepted<FunnelMethod>("funnelsqp"));
}

namespace {
   // a feasible trial iterate that decreases the objective, but with a nonpositive predicted reduction: there is no
   // sufficient decrease to measure it against. The trial iterate is rejected (the trust region shrinks)
   template <typename Strategy>
   bool zero_margin_descent_step_is_accepted(const std::string& preset) {
      const Options options = preset_options(preset);
      Strategy strategy(options);
      Statistics statistics{};
      statistics.add_column("Status", Statistics::string_width, 0);
      statistics.add_column("Funnel", Statistics::double_width, 2);
      Iterate initial_iterate(1, 1);
      initial_iterate.progress = progress(0., 1.);
      strategy.initialize(statistics, initial_iterate);
      return strategy.is_iterate_acceptable(statistics, progress(0., 1.), progress(0., 0.5), progress(0., -0.1), 1.).accepted;
   }
} // namespace

TEST(FeasiblePointAcceptance, NonpositivePredictedReductionAtFeasiblePointIsRejected) {
   EXPECT_FALSE(zero_margin_descent_step_is_accepted<FletcherFilterMethod>("filtersqp"));
   EXPECT_FALSE(zero_margin_descent_step_is_accepted<WaechterFilterMethod>("ipopt"));
   EXPECT_FALSE(zero_margin_descent_step_is_accepted<FunnelMethod>("funnelsqp"));
}

TEST(FeasiblePointAcceptance, HTypeIteratesAtInfeasiblePointsAreUnaffected) {
   // infeasibility halved, objective increased, zero predicted reduction: h-type, accepted
   const Options options = preset_options("filtersqp");
   FletcherFilterMethod strategy(options);
   Statistics statistics{};
   statistics.add_column("Status", Statistics::string_width, 0);
   Iterate initial_iterate(1, 1);
   initial_iterate.progress = progress(1., 1.);
   strategy.initialize(statistics, initial_iterate);
   EXPECT_TRUE(strategy.is_iterate_acceptable(statistics, progress(1., 1.), progress(0.5, 2.), progress(0., 0.), 1.).accepted);
}
