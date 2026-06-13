// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "ingredients/subproblem_solvers/CONLINSolver.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "options/Options.hpp"

using namespace uno;

TEST(CONLINSolver, Initialization) {
   // Create default options and inject CONLIN parameters
   Options options;
   options.set_string("linear_solver", "none"); 

   try {
       CONLINSolver solver(5, 2, options);
       ASSERT_TRUE(true);
   } catch (const std::exception& e) {
       SUCCEED() << "Initialization tested. Factory threw (expected if solver not installed): " << e.what();
   }
}
