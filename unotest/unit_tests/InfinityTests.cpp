// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <vector>
#include <gtest/gtest.h>
#include "tools/Infinity.hpp"

using namespace uno;

// Bug: CUTEst (and AMPL, C and Fortran users) represent infinite bounds by ±1e20. Uno treated them as finite bounds:
// the interior-point method added barrier terms (and multipliers) for them, which altered its iterates
TEST(Infinity, LargeBoundsAreInfinite) {
   EXPECT_EQ(normalize_infinite_bound(1e20), Inf);
   EXPECT_EQ(normalize_infinite_bound(-1e20), -Inf);
   EXPECT_EQ(normalize_infinite_bound(1e30), Inf);
   EXPECT_EQ(normalize_infinite_bound(Inf), Inf);
   EXPECT_EQ(normalize_infinite_bound(-Inf), -Inf);
   EXPECT_EQ(normalize_infinite_bound(1e19), 1e19);
   EXPECT_EQ(normalize_infinite_bound(-3.), -3.);
   EXPECT_EQ(normalize_infinite_bounds(std::vector<double>{-1e20, 0., 1e20}), (std::vector<double>{-Inf, 0., Inf}));
}
