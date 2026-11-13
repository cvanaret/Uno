// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "optimization/WarmstartInformation.hpp"

using namespace uno;

// Bug: iterate_changed() reset the sparsity flags to false. Uno calls it at the beginning of every outer iteration, which
// cancelled the cold start requested by whole_problem_changed() before the first iteration (Uno::uno_solve) and after
// the switch back from the restoration phase (FeasibilityRestoration): the active-set QP solver (BQPD) never received it
// and warm-started from a default/stale active set.
TEST(WarmstartInformation, IterateChangeKeepsPendingStructuralChanges) {
   WarmstartInformation warmstart_information{};
   warmstart_information.whole_problem_changed();
   warmstart_information.iterate_changed();
   EXPECT_TRUE(warmstart_information.new_iterate);
   EXPECT_TRUE(warmstart_information.constraint_bounds_changed);
   EXPECT_TRUE(warmstart_information.trust_region_changed);
   EXPECT_TRUE(warmstart_information.hessian_sparsity_changed);
   EXPECT_TRUE(warmstart_information.jacobian_sparsity_changed);
}

TEST(WarmstartInformation, IterateChangeDoesNotChangeTheStructure) {
   WarmstartInformation warmstart_information{};
   warmstart_information.no_changes();
   warmstart_information.iterate_changed();
   EXPECT_TRUE(warmstart_information.new_iterate);
   EXPECT_TRUE(warmstart_information.constraint_bounds_changed);
   EXPECT_TRUE(warmstart_information.trust_region_changed);
   EXPECT_FALSE(warmstart_information.hessian_sparsity_changed);
   EXPECT_FALSE(warmstart_information.jacobian_sparsity_changed);
}