// Copyright (c) 2018-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ingredients/subproblem_solvers/MUMPS/MUMPSSolver.hpp"
#include "LinearSolverTests.hpp"

namespace uno {
   INSTANTIATE_TYPED_TEST_SUITE_P(MUMPS, LinearSolverTest, MUMPSSolver);
}