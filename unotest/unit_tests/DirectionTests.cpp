// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <limits>
#include <sstream>
#include <gtest/gtest.h>
#include "optimization/Direction.hpp"

using namespace uno;

// Bug: Direction::status_to_string() had no case for SubproblemStatus::ERROR (returned by BQPD and HiGHS upon failure)
// and printed "unknown status, something went wrong" for it
TEST(Direction, PrintsSubproblemSolverError) {
   Direction direction(1, 0);
   direction.status = SubproblemStatus::ERROR;
   std::ostringstream stream;
   stream << direction;
   EXPECT_NE(stream.str().find("status: subproblem solver error"), std::string::npos) << stream.str();
   EXPECT_EQ(stream.str().find("unknown status"), std::string::npos) << stream.str();
}
// A subproblem solver may return a solution with non-finite components (e.g. BQPD on CUTEst NET2, whose Hessian is
// not defined at some iterates): such a direction must be detected and discarded
TEST(Direction, DetectsNonFiniteComponents) {
   Direction direction(2, 1);
   EXPECT_TRUE(direction.is_finite());
   direction.multipliers.constraints[0] = std::numeric_limits<double>::quiet_NaN();
   EXPECT_FALSE(direction.is_finite());
   direction.multipliers.constraints[0] = 0.;
   direction.primals[1] = -std::numeric_limits<double>::infinity();
   EXPECT_FALSE(direction.is_finite());
}
