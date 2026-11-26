// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// The soft restoration phase of the line search (IPOPT, Wächter & Biegler 2006, Section 3.3) is part of the ipopt
// preset. It is disabled by default: it requires an inequality handling method that provides the primal-dual error
// (interior-point methods)

#include <gtest/gtest.h>
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"

using namespace uno;

TEST(SoftRestorationPreset, DisabledByDefault) {
   Options options;
   DefaultOptions::load(options);
   EXPECT_FALSE(options.get_bool("LS_soft_restoration"));
   // IPOPT's soft_resto_pderror_reduction_factor and max_soft_resto_iters
   EXPECT_DOUBLE_EQ(options.get_double("LS_soft_restoration_error_reduction_factor"), 1. - 1e-4);
   EXPECT_EQ(options.get_unsigned_int("LS_soft_restoration_max_iterations"), 10u);
}

TEST(SoftRestorationPreset, EnabledInIpoptPreset) {
   Options options;
   DefaultOptions::load(options);
   Presets::set(options, "ipopt");
   EXPECT_TRUE(options.get_bool("LS_soft_restoration"));
}
