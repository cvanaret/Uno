// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>
#include <gtest/gtest.h>
#include "ingredients/inequality_handling_methods/interior_point_methods/barrier_problems/QualityFunction.hpp"

using namespace uno::quality_function;

namespace {
   // minimum of q on a dense logarithmic grid of [sigma_lo, sigma_up]
   double grid_minimum(const QualityFunction& q, double sigma_lo, double sigma_up, size_t number_points = 20000) {
      double minimum = q(sigma_lo);
      for (size_t k = 0; k <= number_points; ++k) {
         const double sigma = sigma_lo * std::pow(sigma_up / sigma_lo, static_cast<double>(k) / static_cast<double>(number_points));
         minimum = std::min(minimum, q(sigma));
      }
      return minimum;
   }

   // random interior pairs: lower (c, e > 0) and upper (c, e < 0) bounds, with directions of mixed signs
   std::vector<ComplementarityPair> random_pairs(std::mt19937& generator, size_t number_pairs) {
      std::normal_distribution<double> normal(0., 1.);
      std::vector<ComplementarityPair> pairs;
      for (size_t k = 0; k < number_pairs; ++k) {
         const double sign = (k % 2 == 0) ? 1. : -1.;
         const double c = sign * std::exp(2. * normal(generator)), e = sign * std::exp(2. * normal(generator));
         pairs.push_back({c, 2. * c * normal(generator), c * normal(generator), e, 2. * e * normal(generator), e * normal(generator)});
      }
      return pairs;
   }
}

TEST(QualityFunction, PolynomialRoots) {
   // (x - 1)(x - 2)(x - 3) = x^3 - 6x^2 + 11x - 6
   std::vector<double> roots;
   roots_in_interval({-6., 11., -6., 1.}, 0., 4., roots);
   std::sort(roots.begin(), roots.end());
   ASSERT_EQ(roots.size(), 3);
   EXPECT_NEAR(roots[0], 1., 1e-12);
   EXPECT_NEAR(roots[1], 2., 1e-12);
   EXPECT_NEAR(roots[2], 3., 1e-12);
   // roots outside the open interval are discarded
   roots.clear();
   roots_in_interval({-6., 11., -6., 1.}, 1.5, 2.5, roots);
   ASSERT_EQ(roots.size(), 1);
   EXPECT_NEAR(roots[0], 2., 1e-12);
   // no real root
   roots.clear();
   roots_in_interval({1., 0., 1.}, -10., 10., roots);
   EXPECT_TRUE(roots.empty());
}

TEST(QualityFunction, UpperEnvelope) {
   const std::vector<Line> lines{{3., -1.}, {0., 1.}, {1.5, 0.}, {-5., 0.5}};
   const std::vector<Piece> pieces = upper_envelope(lines, 0., 10.);
   ASSERT_FALSE(pieces.empty());
   EXPECT_EQ(pieces.front().start, 0.);
   EXPECT_EQ(pieces.back().end, 10.);
   for (size_t index = 0; index + 1 < pieces.size(); ++index) {
      EXPECT_EQ(pieces[index].end, pieces[index + 1].start);
   }
   // each piece is the maximum of the lines and of the constant line y = 1
   for (const Piece& piece: pieces) {
      const double sigma = 0.5 * (piece.start + piece.end);
      double maximum = 1.;
      for (const Line& line: lines) {
         maximum = std::max(maximum, line.intercept + line.slope * sigma);
      }
      EXPECT_NEAR(piece.line.intercept + piece.line.slope * sigma, maximum, 1e-14);
   }
}

TEST(QualityFunction, NoPairs) {
   const QualityFunction q({}, 0.99, 1., 1., 1.);
   EXPECT_EQ(q(0.5), 0.);
   const Minimizer minimizer = q.minimize(1e-6, 1e2);
   EXPECT_EQ(minimizer.quality, 0.);
}

TEST(QualityFunction, DegenerateInterval) {
   std::mt19937 generator(1);
   const QualityFunction q(random_pairs(generator, 10), 0.99, 1., 1., 0.1);
   const Minimizer minimizer = q.minimize(0.3, 0.3);
   EXPECT_EQ(minimizer.sigma, 0.3);
   EXPECT_EQ(minimizer.quality, q(0.3));
}

TEST(QualityFunction, MinimizerReturnsItsQuality) {
   std::mt19937 generator(2);
   const QualityFunction q(random_pairs(generator, 30), 0.99, 0.5, 2., 1. / 30.);
   const Minimizer minimizer = q.minimize(1e-6, 1e2);
   EXPECT_GE(minimizer.sigma, 1e-6);
   EXPECT_LE(minimizer.sigma, 1e2);
   EXPECT_EQ(minimizer.quality, q(minimizer.sigma));
}

// the exact minimizer is never worse than a dense grid (up to rounding)
TEST(QualityFunction, GlobalMinimumOnRandomInstances) {
   std::mt19937 generator(3);
   std::normal_distribution<double> normal(0., 1.);
   for (size_t instance = 0; instance < 200; ++instance) {
      const size_t number_pairs = 1 + instance % 40;
      const QualityFunction q(random_pairs(generator, number_pairs), 0.99, std::exp(4. * normal(generator)),
         std::exp(4. * normal(generator)), 1. / static_cast<double>(number_pairs));
      const double sigma_lo = 1e-6, sigma_up = 1e2;
      const Minimizer minimizer = q.minimize(sigma_lo, sigma_up);
      const double reference = grid_minimum(q, sigma_lo, sigma_up, 2000);
      EXPECT_LE(minimizer.quality, reference * (1. + 1e-10) + 1e-300) << "instance " << instance;
   }
}

// minimum at a kink where the dual steplength reaches 1, with a huge dual weight (as on CUTEst's TRO41X9): q is very steep
// on one side of the breakpoint, and the rounding of the computed breakpoint lands on that side. Without the points just
// inside the intervals, the minimizer returns 3.7e-5 instead of 1.1e-5
TEST(QualityFunction, KinkWithLargeDualWeight) {
   const ComplementarityPair pair{0.20639983192803246, 0.5804626019152096, 0.20472020428666016, 0.1866973147339866,
      -6.6406735472441989, 1.2932893764119637};
   const double tau = 0.99;
   const QualityFunction q({pair}, tau, 1e24, 0., 1.);
   const double sigma_lo = 1e-6, sigma_up = 1e2;
   const Minimizer minimizer = q.minimize(sigma_lo, sigma_up);
   // the dual steplength reaches 1 at sigma* = (-tau e - g) / h. q is minimal just to its right, and too steep for a grid
   const double kink = (-tau * pair.dual - pair.affine_dual) / pair.centering_dual;
   ASSERT_LT(sigma_lo, kink);
   ASSERT_LT(kink, sigma_up);
   // the candidates just inside the intervals are at a relative distance 1e-10 of the interval width
   EXPECT_LE(minimizer.quality, q(kink * (1. + 1e-9)) * (1. + 1e-4));
   EXPECT_LE(minimizer.quality, grid_minimum(q, sigma_lo, sigma_up) * (1. + 1e-10));
   const auto [primal_step_length, dual_step_length] = q.step_lengths(minimizer.sigma);
   EXPECT_GT(primal_step_length, 0.);
   EXPECT_NEAR(dual_step_length, 1., 1e-8);
}

// many affine pieces (the lines are tangent to a parabola): the cost per interval is O(1)
TEST(QualityFunction, ManyIntervals) {
   const double tau = 0.99;
   const size_t number_pairs = 20000;
   std::vector<ComplementarityPair> pairs;
   for (size_t k = 0; k < number_pairs; ++k) {
      // reciprocal steplength lines 2 t sigma - t^2 + 5 (tangent to sigma^2 + 5, hence all on the envelope)
      const double t = 1. + 9. * static_cast<double>(k) / static_cast<double>(number_pairs);
      const double intercept = 5. - t * t, slope = 2. * t;
      pairs.push_back({1., -tau * intercept, -tau * slope, 1., -tau * intercept, -0.7 * tau * slope});
   }
   const QualityFunction q(pairs, tau, 1., 1., 1. / static_cast<double>(number_pairs));
   const Minimizer minimizer = q.minimize(0.01, 20.);
   EXPECT_GT(minimizer.number_intervals, number_pairs);
   EXPECT_LE(minimizer.quality, grid_minimum(q, 0.01, 20., 500) * (1. + 1e-10));
}
