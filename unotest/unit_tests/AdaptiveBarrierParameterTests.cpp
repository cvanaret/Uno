// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "ingredients/inequality_handling_methods/interior_point_methods/AdaptiveBarrierParameter.hpp"
#include "options/Options.hpp"

using namespace uno;

namespace {
   constexpr double mu_min = 1e-11;
   constexpr double mu_max_factor = 1e3;
   constexpr double average_complementarity = 1.;

   AdaptiveBarrierParameter make_adaptive_barrier_parameter() {
      Options options;
      options.set_double("adaptive_mu_sigma_min", 1e-6);
      options.set_double("adaptive_mu_sigma_max", 1e2);
      options.set_double("adaptive_mu_min", mu_min);
      options.set_double("adaptive_mu_max_factor", mu_max_factor);
      options.set_double("adaptive_mu_monotone_init_factor", 0.8);
      options.set_double("adaptive_mu_filter_margin_factor", 1e-5);
      options.set_double("adaptive_mu_filter_max_margin", 1.);
      return AdaptiveBarrierParameter(options);
   }
}

TEST(AdaptiveBarrierParameter, StartsInFreeMode) {
   AdaptiveBarrierParameter adaptive = make_adaptive_barrier_parameter();
   EXPECT_TRUE(adaptive.in_free_mode());
   EXPECT_EQ(adaptive.update_mode(1., 10., 1., average_complementarity), BarrierUpdate::ORACLE);
}

TEST(AdaptiveBarrierParameter, StaysFreeWhileProgressing) {
   AdaptiveBarrierParameter adaptive = make_adaptive_barrier_parameter();
   EXPECT_EQ(adaptive.update_mode(1., 10., 1., average_complementarity), BarrierUpdate::ORACLE);
   // infeasibility decreases
   EXPECT_EQ(adaptive.update_mode(0.5, 11., 1., average_complementarity), BarrierUpdate::ORACLE);
   // objective decreases
   EXPECT_EQ(adaptive.update_mode(0.6, 9., 1., average_complementarity), BarrierUpdate::ORACLE);
   EXPECT_TRUE(adaptive.in_free_mode());
}

TEST(AdaptiveBarrierParameter, SwitchesToFixedModeThenMonotone) {
   AdaptiveBarrierParameter adaptive = make_adaptive_barrier_parameter();
   EXPECT_EQ(adaptive.update_mode(1., 10., 1., average_complementarity), BarrierUpdate::ORACLE);
   // dominated by the filter entry (1, 10)
   EXPECT_EQ(adaptive.update_mode(1.5, 12., 1., average_complementarity), BarrierUpdate::RESET_FIXED);
   EXPECT_FALSE(adaptive.in_free_mode());
   EXPECT_EQ(adaptive.update_mode(1.5, 12., 1., average_complementarity), BarrierUpdate::MONOTONE);
   // progress again: back to the free mode
   EXPECT_EQ(adaptive.update_mode(0.1, 12., 1., average_complementarity), BarrierUpdate::ORACLE);
   EXPECT_TRUE(adaptive.in_free_mode());
}

TEST(AdaptiveBarrierParameter, FilterMargin) {
   AdaptiveBarrierParameter adaptive = make_adaptive_barrier_parameter();
   EXPECT_EQ(adaptive.update_mode(1., 10., 1., average_complementarity), BarrierUpdate::ORACLE);
   // a decrease smaller than the margin 1e-5 * min(1, error) is not sufficient progress
   EXPECT_EQ(adaptive.update_mode(1. - 1e-6, 10. - 1e-6, 1., average_complementarity), BarrierUpdate::RESET_FIXED);
}

TEST(AdaptiveBarrierParameter, LeaveFreeModeAndReset) {
   AdaptiveBarrierParameter adaptive = make_adaptive_barrier_parameter();
   EXPECT_EQ(adaptive.update_mode(1., 10., 1., average_complementarity), BarrierUpdate::ORACLE);
   adaptive.leave_free_mode();
   EXPECT_FALSE(adaptive.in_free_mode());
   EXPECT_EQ(adaptive.update_mode(2., 20., 1., average_complementarity), BarrierUpdate::MONOTONE);
   // after the restoration phase: free mode with an empty filter (any point makes progress)
   adaptive.reset();
   EXPECT_TRUE(adaptive.in_free_mode());
   EXPECT_EQ(adaptive.update_mode(2., 20., 1., average_complementarity), BarrierUpdate::ORACLE);
}

TEST(AdaptiveBarrierParameter, BarrierParameterBounds) {
   AdaptiveBarrierParameter adaptive = make_adaptive_barrier_parameter();
   // the upper bound mu_max_factor * average complementarity is set at the first iterate
   (void) adaptive.update_mode(1., 10., 1., average_complementarity);
   const double mu_max = mu_max_factor * average_complementarity;
   EXPECT_EQ(adaptive.clamp(1e-20), mu_min);
   EXPECT_EQ(adaptive.clamp(1e10), mu_max);
   EXPECT_EQ(adaptive.clamp(0.5), 0.5);
   EXPECT_EQ(adaptive.fixed_mode_barrier_parameter(0.5), 0.8 * 0.5);
   // the average complementarity at later iterates does not change mu_max
   (void) adaptive.update_mode(0.5, 10., 1., 1e-3);
   EXPECT_EQ(adaptive.clamp(1e10), mu_max);
}

TEST(AdaptiveBarrierParameter, SigmaInterval) {
   AdaptiveBarrierParameter adaptive = make_adaptive_barrier_parameter();
   (void) adaptive.update_mode(1., 10., 1., average_complementarity);
   const auto [sigma_lo, sigma_up] = adaptive.sigma_interval(1.);
   EXPECT_EQ(sigma_lo, 1e-6);
   EXPECT_EQ(sigma_up, 1e2);
   // sigma * average_complementarity is kept within [mu_min, mu_max]
   const auto [tiny_lo, tiny_up] = adaptive.sigma_interval(1e-9);
   EXPECT_DOUBLE_EQ(tiny_lo, mu_min / 1e-9);
   EXPECT_EQ(tiny_up, 1e2);
   const auto [large_lo, large_up] = adaptive.sigma_interval(1e6);
   EXPECT_EQ(large_lo, 1e-6);
   EXPECT_DOUBLE_EQ(large_up, 1e3 / 1e6);
   // empty interval: collapsed to its lower bound
   const auto [empty_lo, empty_up] = adaptive.sigma_interval(1e12);
   EXPECT_EQ(empty_lo, empty_up);
}
