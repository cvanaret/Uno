// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <stdexcept>
#include <gtest/gtest.h>
#include "ingredients/inequality_handling_methods/interior_point_methods/BarrierParameterUpdateStrategy.hpp"
#include "options/Options.hpp"

using namespace uno;

namespace {
   // only update_barrier_parameter() uses the barrier problem, so the other member functions can be tested with a dummy
   struct DummyBarrierProblem { };
   using Strategy = BarrierParameterUpdateStrategy<DummyBarrierProblem>;

   constexpr double initial_parameter = 0.1;
   constexpr double dual_tolerance = 1e-8;
   constexpr double k_mu = 0.2;
   constexpr double theta_mu = 1.5;
   constexpr double k_epsilon = 10.;
   constexpr double update_fraction = 10.;
   constexpr double minimal_parameter = dual_tolerance / update_fraction;

   Options make_barrier_options() {
      Options options;
      options.set_double("barrier_initial_parameter", initial_parameter);
      options.set_double("dual_tolerance", dual_tolerance);
      options.set_double("barrier_k_mu", k_mu);
      options.set_double("barrier_theta_mu", theta_mu);
      options.set_double("barrier_k_epsilon", k_epsilon);
      options.set_double("barrier_update_fraction", update_fraction);
      return options;
   }
}

TEST(BarrierParameterUpdateStrategy, InitialParameter) {
   const Strategy strategy(make_barrier_options());
   EXPECT_EQ(strategy.get_barrier_parameter(), initial_parameter);
}

TEST(BarrierParameterUpdateStrategy, SetRejectsNonPositiveParameter) {
   Strategy strategy(make_barrier_options());
   EXPECT_THROW(strategy.set_barrier_parameter(0.), std::runtime_error);
   EXPECT_THROW(strategy.set_barrier_parameter(-1.), std::runtime_error);
   EXPECT_EQ(strategy.get_barrier_parameter(), initial_parameter);
   constexpr double new_barrier_parameter = 1e-3;
   strategy.set_barrier_parameter(new_barrier_parameter);
   EXPECT_EQ(strategy.get_barrier_parameter(), new_barrier_parameter);
}

// mu > k_mu^(1/(theta_mu - 1)) = 0.04: linear decrease k_mu * mu
TEST(BarrierParameterUpdateStrategy, ForcedDecreaseLinearRegime) {
   Strategy strategy(make_barrier_options()); // mu = 0.1
   EXPECT_TRUE(strategy.force_barrier_parameter_decrease());
   EXPECT_DOUBLE_EQ(strategy.get_barrier_parameter(), k_mu * initial_parameter);
}

// mu < k_mu^(1/(theta_mu - 1)) = 0.04: superlinear decrease mu^theta_mu
TEST(BarrierParameterUpdateStrategy, ForcedDecreaseSuperlinearRegime) {
   Strategy strategy(make_barrier_options());
   strategy.set_barrier_parameter(1e-2);
   EXPECT_TRUE(strategy.force_barrier_parameter_decrease());
   EXPECT_DOUBLE_EQ(strategy.get_barrier_parameter(), std::pow(1e-2, theta_mu));
}

// stress test: repeated forced decreases are strictly monotone, terminate, and saturate at dual_tolerance/update_fraction
TEST(BarrierParameterUpdateStrategy, ForcedDecreaseSaturatesAtLowerBound) {
   Strategy strategy(make_barrier_options());
   strategy.set_barrier_parameter(1e6);
   size_t number_decreases = 0;
   double previous_parameter = strategy.get_barrier_parameter();
   while (strategy.force_barrier_parameter_decrease()) {
      ASSERT_LT(strategy.get_barrier_parameter(), previous_parameter);
      ASSERT_GE(strategy.get_barrier_parameter(), minimal_parameter);
      previous_parameter = strategy.get_barrier_parameter();
      ASSERT_LT(++number_decreases, 100u);
   }
   EXPECT_EQ(strategy.get_barrier_parameter(), minimal_parameter);
   EXPECT_FALSE(strategy.force_barrier_parameter_decrease());
   EXPECT_EQ(strategy.get_barrier_parameter(), minimal_parameter);
}

// a parameter set below the lower bound is pulled back up by a forced "decrease", which reports no decrease
TEST(BarrierParameterUpdateStrategy, ForcedDecreaseBelowLowerBound) {
   Strategy strategy(make_barrier_options());
   strategy.set_barrier_parameter(minimal_parameter / 100.);
   EXPECT_FALSE(strategy.force_barrier_parameter_decrease());
   EXPECT_EQ(strategy.get_barrier_parameter(), minimal_parameter);
}