// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <limits>
#include <stdexcept>
#include <gtest/gtest.h>
#include "linear_algebra/Norm.hpp"
#include "linear_algebra/Vector.hpp"
#include "symbolic/Range.hpp"

using namespace uno;

TEST(Norm, SingleVector) {
   const Vector<double> x{3., -4.};
   EXPECT_DOUBLE_EQ(norm_1(x), 7.);
   EXPECT_DOUBLE_EQ(norm_2_squared(x), 25.);
   EXPECT_DOUBLE_EQ(norm_2(x), 5.);
   EXPECT_DOUBLE_EQ(norm_inf(x), 4.);
}

TEST(Norm, SeveralVectors) {
   const Vector<double> x{3., -4.};
   const Vector<double> y{-12.};
   EXPECT_DOUBLE_EQ(norm_1(x, y), 19.);
   EXPECT_DOUBLE_EQ(norm_2_squared(x, y), 169.);
   EXPECT_DOUBLE_EQ(norm_2(x, y), 13.);
   EXPECT_DOUBLE_EQ(norm_inf(x, y), 12.);
}

TEST(Norm, Dispatch) {
   const Vector<double> x{1., -2., 2.};
   EXPECT_DOUBLE_EQ(norm(Norm::L1, x), 5.);
   EXPECT_DOUBLE_EQ(norm(Norm::L2, x), 3.);
   EXPECT_DOUBLE_EQ(norm(Norm::L2_SQUARED, x), 9.);
   EXPECT_DOUBLE_EQ(norm(Norm::INF, x), 2.);
}

TEST(Norm, EmptyVectorHasZeroNorm) {
   const Vector<double> x(0);
   EXPECT_EQ(norm_1(x), 0.);
   EXPECT_EQ(norm_2(x), 0.);
   EXPECT_EQ(norm_inf(x), 0.);
}

TEST(Norm, FromString) {
   EXPECT_EQ(norm_from_string("L1"), Norm::L1);
   EXPECT_EQ(norm_from_string("L2"), Norm::L2);
   EXPECT_EQ(norm_from_string("L2_squared"), Norm::L2_SQUARED);
   EXPECT_EQ(norm_from_string("INF"), Norm::INF);
   EXPECT_THROW(norm_from_string("l1"), std::invalid_argument); // case sensitive
   EXPECT_THROW(norm_from_string(""), std::invalid_argument);
}

TEST(Norm, InfinityPropagates) {
   const Vector<double> x{1., -INFINITY};
   EXPECT_EQ(norm_1(x), INFINITY);
   EXPECT_EQ(norm_2(x), INFINITY);
   EXPECT_EQ(norm_inf(x), INFINITY);
}

TEST(Norm, NaNPropagates) {
   const Vector<double> x{1., std::numeric_limits<double>::quiet_NaN()};
   EXPECT_TRUE(std::isnan(norm_1(x)));
   EXPECT_TRUE(std::isnan(norm_2(x)));
   EXPECT_TRUE(std::isnan(norm_inf(x)));
}

// stress test: equivalence inequalities ||x||_inf <= ||x||_2 <= ||x||_1 <= n ||x||_inf on many vectors
TEST(Norm, EquivalenceInequalities) {
   for (size_t size: {1u, 2u, 10u, 1000u}) {
      Vector<double> x(size);
      for (size_t index: Range(size)) {
         x[index] = std::sin(1. + 37. * static_cast<double>(index)) * std::pow(10., static_cast<double>(index % 7) - 3.);
      }
      const double norm_1_x = norm_1(x);
      const double norm_2_x = norm_2(x);
      const double norm_inf_x = norm_inf(x);
      const double tolerance = 1e-12 * norm_1_x;
      EXPECT_LE(norm_inf_x, norm_2_x + tolerance);
      EXPECT_LE(norm_2_x, norm_1_x + tolerance);
      EXPECT_LE(norm_1_x, static_cast<double>(size) * norm_inf_x + tolerance);
   }
}