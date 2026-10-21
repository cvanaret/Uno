// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

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