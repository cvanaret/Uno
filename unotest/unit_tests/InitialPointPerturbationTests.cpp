// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <vector>
#include <gtest/gtest.h>
#include "optimization/InitialPointPerturbation.hpp"
#include "tools/Infinity.hpp"

using namespace uno;

TEST(InitialPointPerturbation, IsDeterministic) {
   const std::vector<double> lower_bounds(4, -Inf), upper_bounds(4, Inf);
   Vector<double> x1{0., 1., -3., 100.};
   Vector<double> x2{0., 1., -3., 100.};
   perturb_point(x1, lower_bounds, upper_bounds, 1e-6);
   perturb_point(x2, lower_bounds, upper_bounds, 1e-6);
   for (size_t index = 0; index < 4; ++index) {
      EXPECT_EQ(x1[index], x2[index]);
   }
}

TEST(InitialPointPerturbation, MovesEachFreeVariableByTheRelativeMagnitude) {
   const double magnitude = 1e-4;
   const std::vector<double> lower_bounds(3, -Inf), upper_bounds(3, Inf);
   const Vector<double> original{0., 0.5, -200.};
   Vector<double> x = original;
   perturb_point(x, lower_bounds, upper_bounds, magnitude);
   for (size_t index = 0; index < 3; ++index) {
      const double distance = std::abs(x[index] - original[index]);
      const double scale = magnitude * std::max(1., std::abs(original[index]));
      EXPECT_GE(distance, scale);
      EXPECT_LT(distance, 2. * scale);
   }
}

// coincident coordinates (e.g. a distance between two points evaluated at 0) are separated
TEST(InitialPointPerturbation, SeparatesCoincidentCoordinates) {
   const size_t n = 50;
   const std::vector<double> lower_bounds(n, -Inf), upper_bounds(n, Inf);
   Vector<double> x(n, 0.);
   perturb_point(x, lower_bounds, upper_bounds, 1e-6);
   for (size_t i = 0; i < n; ++i) {
      for (size_t j = i + 1; j < n; ++j) {
         EXPECT_NE(x[i], x[j]) << "indices " << i << " and " << j;
      }
   }
}

TEST(InitialPointPerturbation, StaysWithinBoundsAndKeepsFixedVariables) {
   const std::vector<double> lower_bounds{0., -1., 2., -Inf, 3.};
   const std::vector<double> upper_bounds{1., 0., 2., 5., Inf};
   // at the lower bound, at the upper bound, fixed, at the upper bound, at the lower bound
   Vector<double> x{0., 0., 2., 5., 3.};
   perturb_point(x, lower_bounds, upper_bounds, 1e-3);
   for (size_t index = 0; index < x.size(); ++index) {
      EXPECT_GE(x[index], lower_bounds[index]);
      EXPECT_LE(x[index], upper_bounds[index]);
   }
   EXPECT_EQ(x[2], 2.); // fixed variable untouched
   // a variable at one of its bounds is moved into the interior
   EXPECT_GT(x[0], 0.);
   EXPECT_LT(x[1], 0.);
   EXPECT_LT(x[3], 5.);
   EXPECT_GT(x[4], 3.);
}
