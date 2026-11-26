// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Detection of an unbounded objective along feasible iterates: absolute threshold, or decrease by many orders of
// magnitude relative to the first feasible iterate, sustained for several consecutive feasible iterates

#include <cstddef>
#include <gtest/gtest.h>
#include "ingredients/constraint_relaxation_strategies/UnboundedObjectiveDetector.hpp"
#include "tools/Infinity.hpp"

using namespace uno;

namespace {
   // absolute threshold -1e20, relative decrease 1e10, 3 consecutive iterates
   UnboundedObjectiveDetector make_detector() {
      return UnboundedObjectiveDetector(-1e20, 1e10, 3);
   }
} // namespace

// below the absolute threshold: unbounded at once
TEST(UnboundedObjective, AbsoluteThreshold) {
   UnboundedObjectiveDetector detector = make_detector();
   EXPECT_TRUE(detector.is_unbounded(-1e21));
}

// CUTE coshfun (preset ipopt): the iterates follow the feasible ray f = 5 - 2t^2 from a first feasible objective of
// order 1, and the objective decreases by about 1e8 per iteration: -1e20 is out of reach, -1e10 is reached within a few
// hundred iterations
TEST(UnboundedObjective, RelativeDecreaseSustainedOverConsecutiveIterates) {
   UnboundedObjectiveDetector detector = make_detector();
   EXPECT_FALSE(detector.is_unbounded(0.2)); // first feasible iterate: reference
   EXPECT_FALSE(detector.is_unbounded(-1e9));
   EXPECT_FALSE(detector.is_unbounded(-2e10)); // 1 of 3
   EXPECT_FALSE(detector.is_unbounded(-3e10)); // 2 of 3
   EXPECT_TRUE(detector.is_unbounded(-4e10));  // 3 of 3
}

// the streak counts consecutive feasible iterates: an infeasible iterate or an iterate above the threshold restarts it
TEST(UnboundedObjective, StreakIsInterrupted) {
   UnboundedObjectiveDetector detector = make_detector();
   EXPECT_FALSE(detector.is_unbounded(0.));
   EXPECT_FALSE(detector.is_unbounded(-2e10));
   EXPECT_FALSE(detector.is_unbounded(-3e10));
   detector.notify_infeasible_iterate();
   EXPECT_FALSE(detector.is_unbounded(-4e10));
   EXPECT_FALSE(detector.is_unbounded(-5e10));
   EXPECT_FALSE(detector.is_unbounded(-1.)); // above the threshold
   EXPECT_FALSE(detector.is_unbounded(-6e10));
   EXPECT_FALSE(detector.is_unbounded(-7e10));
   EXPECT_TRUE(detector.is_unbounded(-8e10));
}

// the threshold is relative to the magnitude of the reference objective (badly scaled models)
TEST(UnboundedObjective, ThresholdIsRelativeToReferenceObjective) {
   UnboundedObjectiveDetector detector = make_detector();
   EXPECT_FALSE(detector.is_unbounded(-1e5)); // reference: threshold = -1e5 - 1e10 * 1e5 ~ -1e15
   for (size_t iteration = 0; iteration < 10; ++iteration) {
      EXPECT_FALSE(detector.is_unbounded(-1e14));
   }
   EXPECT_FALSE(detector.is_unbounded(-2e15));
   EXPECT_FALSE(detector.is_unbounded(-2e15));
   EXPECT_TRUE(detector.is_unbounded(-2e15));
}

// an infinite relative decrease disables the relative test
TEST(UnboundedObjective, RelativeTestDisabled) {
   UnboundedObjectiveDetector detector(-1e20, Inf, 3);
   EXPECT_FALSE(detector.is_unbounded(0.));
   for (size_t iteration = 0; iteration < 10; ++iteration) {
      EXPECT_FALSE(detector.is_unbounded(-1e19));
   }
   EXPECT_TRUE(detector.is_unbounded(-1e21));
}
