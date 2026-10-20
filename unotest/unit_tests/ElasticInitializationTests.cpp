// Copyright (c) 2018-2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include <cmath>
#include "ingredients/inequality_handling_methods/interior_point_methods/ElasticInitialization.hpp"

using namespace uno;

namespace {
   void check_pair(double r, double epsilon) {
      const auto [a, b] = solve_elastic_pair(r, epsilon);
      ASSERT_GT(a, 0.);
      ASSERT_GT(b, 0.);
      ASSERT_GE(std::min(a, b), epsilon / 2.);
      // a - b = r (relative to the magnitude of the terms)
      EXPECT_NEAR(a - b, r, 1e-14 * std::max(1., std::max(a, b)));
      // 1/a + 1/b = 2/epsilon
      EXPECT_NEAR(epsilon / a + epsilon / b, 2., 1e-12);
   }
}

TEST(ElasticInitialization, Feasible) {
   const auto [a, b] = solve_elastic_pair(0., 1e-3);
   EXPECT_DOUBLE_EQ(a, 1e-3);
   EXPECT_DOUBLE_EQ(b, 1e-3);
}

TEST(ElasticInitialization, Identities) {
   for (const double epsilon: {1e-12, 1e-8, 1e-3, 1., 10.}) {
      for (const double r: {-1e8, -1e2, -1., -1e-6, 0., 1e-6, 1., 1e2, 1e8}) {
         check_pair(r, epsilon);
      }
   }
}

TEST(ElasticInitialization, NoCancellation) {
   // |r| >> epsilon: the naive formula (epsilon - r + sqrt(r^2 + epsilon^2))/2 loses most digits
   const double epsilon = 1e-9;
   const double r = 1e4;
   const auto [a, b] = solve_elastic_pair(r, epsilon);
   EXPECT_DOUBLE_EQ(a, r + b);
   // b = epsilon/2 + epsilon^2/(4r) + O(epsilon^4/r^3)
   EXPECT_DOUBLE_EQ(b, epsilon / 2. + epsilon * epsilon / (4. * r));
   const double naive_b = (epsilon - r + std::sqrt(r * r + epsilon * epsilon)) / 2.;
   EXPECT_GT(std::abs(naive_b - b) / b, 1e-6);
   const auto [c, d] = solve_elastic_pair(-r, epsilon);
   EXPECT_DOUBLE_EQ(c, b);
   EXPECT_DOUBLE_EQ(d, a);
}
