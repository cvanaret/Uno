// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <limits>
#include <gtest/gtest.h>
#include "ingredients/inequality_handling_methods/interior_point_methods/barrier_problems/PrimalDualInteriorPointProblem.hpp"

using namespace uno;

namespace {
   const double eps = std::numeric_limits<double>::epsilon();
   const double slack_move = std::pow(eps, 0.75);

   double safe_slack(double slack, double multiplier, double mu, double bound) {
      return PrimalDualInteriorPointProblem::compute_safe_slack(slack, multiplier, mu, bound, slack_move);
   }
}

TEST(SafeSlack, LargeEnoughSlackIsUnchanged) {
   EXPECT_EQ(safe_slack(1e-3, 1., 0.1, 0.), 1e-3);
   // threshold is eps*min(1, mu)
   EXPECT_EQ(safe_slack(eps, 1., 10., 0.), eps);
   EXPECT_EQ(safe_slack(eps * 1e-2, 1., 1e-2, 0.), eps * 1e-2);
}

TEST(SafeSlack, CenteredSlackWhenBelowCap) {
   // mu/z = 1e-14 lies between s_min and the cap slack_move*max(1, |b|) ~ 1.8e-12
   const double mu = 1e-2, z = 1e12;
   EXPECT_DOUBLE_EQ(safe_slack(0., z, mu, 0.), mu / z);
   // the sign of the multiplier is irrelevant (upper bounds)
   EXPECT_DOUBLE_EQ(safe_slack(0., -z, mu, 0.), mu / z);
}

TEST(SafeSlack, CappedBySlackMove) {
   // mu/z = 1e-1 >> cap
   EXPECT_DOUBLE_EQ(safe_slack(0., 1., 0.1, 0.), slack_move);
   // the cap scales with max(1, |bound|)
   EXPECT_DOUBLE_EQ(safe_slack(0., 1., 0.1, -1e4), slack_move * 1e4);
   // zero multiplier: centered slack is infinite, the cap applies
   EXPECT_DOUBLE_EQ(safe_slack(0., 0., 0.1, 0.), slack_move);
}

TEST(SafeSlack, FlooredByMinimumSlack) {
   // mu/z << s_min = eps*min(1, mu)
   const double mu = 0.5;
   EXPECT_DOUBLE_EQ(safe_slack(0., 1e300, mu, 0.), eps * mu);
}

TEST(SafeSlack, NegativeSlackBecomesPositive) {
   // iterate outside the bound (e.g. after a phase switch): the result is positive
   EXPECT_GT(safe_slack(-1e-6, 1., 0.1, 0.), 0.);
   EXPECT_DOUBLE_EQ(safe_slack(-1e-6, 1., 0.1, 0.), slack_move);
   EXPECT_GE(safe_slack(-1., 1e300, 0.1, 0.), eps * 0.1);
}

TEST(SafeSlack, PartialSlackIsKeptInCap) {
   // 0 < slack < s_min: cap = slack + slack_move*max(1, |b|)
   const double mu = 1., slack = 0.5 * eps;
   EXPECT_DOUBLE_EQ(safe_slack(slack, 1., mu, 0.), slack + slack_move);
}
