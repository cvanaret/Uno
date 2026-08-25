// Copyright (c) 2018-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ingredients/subproblem_solvers/MA57/MA57Solver.hpp"
#include "LinearSolverTests.hpp"

namespace uno {
   INSTANTIATE_TYPED_TEST_SUITE_P(MA57, LinearSolverTest, MA57Solver);
}