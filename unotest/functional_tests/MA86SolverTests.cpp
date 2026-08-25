// Copyright (c) 2026 Alexis Montoison and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ingredients/subproblem_solvers/MA86/MA86Solver.hpp"
#include "LinearSolverTests.hpp"

namespace uno {
   INSTANTIATE_TYPED_TEST_SUITE_P(MA86, LinearSolverTest, MA86Solver);
}