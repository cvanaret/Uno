// Copyright (c) 2018-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include <array>
#include "ingredients/subproblem_solvers/MA27/MA27Solver.hpp"
#include "linear_algebra/Vector.hpp"

using namespace uno;

TEST(MA27Solver, SystemSize5) {
   MA27Solver solver;
   COOLinearSystem& linear_system = solver.get_coo_linear_system();
   constexpr size_t dimension = 5;
   linear_system.dimension = dimension;
   linear_system.number_nonzeros = 7;
   // indices with Fortran-based indexing
   linear_system.matrix_row_indices = {1, 1, 2, 2, 3, 3, 5};
   linear_system.matrix_column_indices = {1, 2, 3, 5, 3, 4, 5};
   linear_system.matrix_values = {2., 3., 4., 6., 1., 5., 1.};
   linear_system.rhs = {8., 45., 31., 15., 17.};
   solver.initialize_memory();
   solver.do_symbolic_analysis();
   solver.do_numerical_factorization(false);
   Vector<double> result(dimension);
   result.fill(0.);
   solver.solve_indefinite_system(result.data());

   const std::array<double, dimension> reference{1., 2., 3., 4., 5.};
   const double tolerance = 1e-8;
   for (size_t index: Range(dimension)) {
      EXPECT_NEAR(result[index], reference[index], tolerance);
   }
}

TEST(MA27Solver, Inertia) {
   MA27Solver solver;
   COOLinearSystem& linear_system = solver.get_coo_linear_system();
   linear_system.dimension = 5;
   linear_system.number_nonzeros = 7;
   linear_system.matrix_row_indices = {1, 1, 2, 2, 3, 3, 5};
   linear_system.matrix_column_indices = {1, 2, 3, 5, 3, 4, 5};
   linear_system.matrix_values = {2., 3., 4., 6., 1., 5., 1.};
   solver.initialize_memory();
   solver.do_symbolic_analysis();
   solver.do_numerical_factorization(false);
   const auto [number_positive, number_negative, number_zero] = solver.get_inertia();
   ASSERT_EQ(number_positive, 3);
   ASSERT_EQ(number_negative, 2);
   ASSERT_EQ(number_zero, 0);
}

TEST(MA27Solver, SingularMatrix) {
   // comes from hs015 solved with byrd preset
   MA27Solver solver;
   COOLinearSystem& linear_system = solver.get_coo_linear_system();
   linear_system.dimension = 4;
   linear_system.number_nonzeros = 7;
   linear_system.matrix_row_indices = {1, 1, 1, 2, 2, 3, 4};
   linear_system.matrix_column_indices = {1, 1, 2, 2, 2, 3, 4};
   linear_system.matrix_values = {-0.0198, 0.625075, -0.277512, -0.624975, 0.625075, 0., 0.};
   solver.initialize_memory();
   solver.do_symbolic_analysis();
   solver.do_numerical_factorization(false);

   // expected inertia (1, 1, 2)
   ASSERT_TRUE(solver.matrix_is_singular());
}