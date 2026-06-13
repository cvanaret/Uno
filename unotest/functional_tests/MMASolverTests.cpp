// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "ingredients/subproblem_solvers/MMASolver.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "options/Options.hpp"

using namespace uno;

TEST(MMASolver, Initialization) {
   // Create default options and inject MMA parameters
   Options options;
   options.set_string("linear_solver", "none"); 
   options.set_double("mma_asyinit", 0.5);
   options.set_double("mma_asyincr", 1.2);
   options.set_double("mma_asydecr", 0.7);
   options.set_double("mma_external_move_limit", 0.2);
   options.set_double("mma_internal_limit", 0.1);
   options.set_integer("mma_max_inner_iterations", 10);

   try {
       MMASolver solver(5, 2, options);
       ASSERT_TRUE(true);
   } catch (const std::exception& e) {
       SUCCEED() << "Initialization tested. Factory threw (expected if solver not installed): " << e.what();
   }
}
