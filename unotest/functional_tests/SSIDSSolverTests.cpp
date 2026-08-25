// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ingredients/subproblem_solvers/SSIDS/SSIDSSolver.hpp"
#include "LinearSolverTests.hpp"

namespace uno {
   INSTANTIATE_TYPED_TEST_SUITE_P(SSIDS, LinearSolverTest, SSIDSSolver);
}