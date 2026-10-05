// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <random>
#include <stdexcept>
#include <gtest/gtest.h>
#include "ingredients/globalization_strategies/switching_methods/funnel_methods/Funnel.hpp"
#include "options/Options.hpp"

using namespace uno;

namespace {
   constexpr double funnel_beta = 0.5;
   constexpr double funnel_kappa = 0.5;

   Options make_funnel_options(uno_int update_strategy, double beta = funnel_beta, double kappa = funnel_kappa) {
      Options options;
      options.set_double("funnel_beta", beta);
      options.set_double("funnel_kappa", kappa);
      options.set_integer("funnel_update_strategy", update_strategy);
      return options;
   }
}

TEST(Funnel, Acceptability) {
   Funnel funnel(make_funnel_options(1));
   funnel.set_infeasibility_upper_bound(1.);
   EXPECT_TRUE(funnel.acceptable(1.));
   EXPECT_TRUE(funnel.acceptable(0.));
   EXPECT_FALSE(funnel.acceptable(1. + 1e-12));
}

TEST(Funnel, SufficientDecreaseUsesMargin) {
   Funnel funnel(make_funnel_options(1));
   funnel.set_infeasibility_upper_bound(1.);
   EXPECT_TRUE(funnel.sufficient_decrease_condition(funnel_beta));
   EXPECT_FALSE(funnel.sufficient_decrease_condition(0.75)); // acceptable, but not a sufficient decrease
   EXPECT_TRUE(funnel.acceptable(0.75));
}

// strategy 1, trial <= current: width = max(beta * width, kappa * current + (1 - kappa) * trial)
TEST(Funnel, Strategy1ConvexCombinationWins) {
   Funnel funnel(make_funnel_options(1));
   funnel.set_infeasibility_upper_bound(1.);
   funnel.update(/* current_infeasibility = */ 0.9, /* trial_infeasibility = */ 0.5);
   EXPECT_DOUBLE_EQ(funnel.current_width(), 0.7);
}

TEST(Funnel, Strategy1MarginWins) {
   Funnel funnel(make_funnel_options(1));
   funnel.set_infeasibility_upper_bound(1.);
   funnel.update(/* current_infeasibility = */ 0.2, /* trial_infeasibility = */ 0.);
   EXPECT_DOUBLE_EQ(funnel.current_width(), funnel_beta);
}

TEST(Funnel, Strategy1InfeasibilityIncreaseShrinksByMargin) {
   Funnel funnel(make_funnel_options(1));
   funnel.set_infeasibility_upper_bound(1.);
   funnel.update(/* current_infeasibility = */ 0.2, /* trial_infeasibility = */ 0.3);
   EXPECT_DOUBLE_EQ(funnel.current_width(), funnel_beta);
}

TEST(Funnel, Strategy2) {
   Funnel funnel(make_funnel_options(2));
   funnel.set_infeasibility_upper_bound(1.);
   funnel.update(/* current_infeasibility = */ 0.9, /* trial_infeasibility = */ 0.2);
   EXPECT_DOUBLE_EQ(funnel.current_width(), 0.6); // 0.5 * 1 + 0.5 * 0.2
}

TEST(Funnel, Strategy3) {
   Funnel funnel(make_funnel_options(3));
   funnel.set_infeasibility_upper_bound(1.);
   funnel.update(/* current_infeasibility = */ 0.9, /* trial_infeasibility = */ 0.2);
   EXPECT_DOUBLE_EQ(funnel.current_width(), funnel_beta);
}

TEST(Funnel, UnknownStrategyThrowsOnUpdate) {
   Funnel funnel(make_funnel_options(42));
   funnel.set_infeasibility_upper_bound(1.);
   EXPECT_THROW(funnel.update(0.5, 0.1), std::runtime_error);
   EXPECT_DOUBLE_EQ(funnel.current_width(), 1.);
}

TEST(Funnel, RestorationUpdate) {
   Funnel funnel(make_funnel_options(1));
   funnel.set_infeasibility_upper_bound(1.);
   funnel.update_restoration(/* current_infeasibility = */ 0.4);
   EXPECT_DOUBLE_EQ(funnel.current_width(), 0.7); // 0.5 * 1 + 0.5 * 0.4
}

// without an upper bound, the width is infinite and no update strategy can shrink it
TEST(Funnel, WidthStaysInfiniteWithoutUpperBound) {
   for (uno_int strategy: {1, 2, 3}) {
      Funnel funnel(make_funnel_options(strategy));
      funnel.update(1., 0.5);
      EXPECT_EQ(funnel.current_width(), Inf) << "strategy " << strategy;
   }
}

TEST(Funnel, Strategy3ShrinksGeometrically) {
   Funnel funnel(make_funnel_options(3));
   funnel.set_infeasibility_upper_bound(1.);
   double expected_width = 1.;
   for (size_t iteration = 0; iteration < 100; ++iteration) {
      funnel.update(/* current_infeasibility = */ 0., /* trial_infeasibility = */ 0.);
      expected_width *= funnel_beta;
      ASSERT_DOUBLE_EQ(funnel.current_width(), expected_width);
   }
}

// stress test: for strategy 1 with trial <= current <= width, the new width lies in [max(beta*width, trial), width],
// so the funnel never grows and the accepted trial point remains inside it
TEST(Funnel, Strategy1RandomizedInvariants) {
   std::mt19937 generator(12345);
   std::uniform_real_distribution<double> unit(0., 1.);
   for (double beta: {0.1, 0.9, 0.9999}) {
      for (double kappa: {0., 0.5, 1.}) {
         Funnel funnel(make_funnel_options(1, beta, kappa));
         double width = 1e3;
         funnel.set_infeasibility_upper_bound(width);
         for (size_t iteration = 0; iteration < 1000; ++iteration) {
            const double current = width * unit(generator); // <= width
            const double trial = current * unit(generator); // <= current
            funnel.update(/* current_infeasibility = */ current, /* trial_infeasibility = */ trial);
            const double new_width = funnel.current_width();
            ASSERT_LE(new_width, width);
            ASSERT_GE(new_width, beta * width * (1. - 1e-15));
            ASSERT_TRUE(funnel.acceptable(trial));
            width = new_width;
         }
      }
   }
}