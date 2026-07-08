// Copyright (c) 2026 Alexis Montoison and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include <array>
#include "ingredients/subproblem_solvers/MA86/MA86Solver.hpp"
#include "linear_algebra/Indexing.hpp"
#include "linear_algebra/Vector.hpp"

using namespace uno;

// MA86 expects the lower triangle of the matrix (row >= column) in C/0-based indexing. These are the
// same symmetric systems used in the MA57 tests, transposed to the lower triangle and 0-indexed.
// The *FortranIndexing tests below feed the same systems in 1-based indexing to check that
// MA86Solver builds a 1-based CSC and drives MA86/MC68 through their C interface (f_arrays=1).

TEST(MA86Solver, SystemSize5) {
   MA86Solver solver;
   COOLinearSystem& linear_system = solver.get_coo_linear_system();
   constexpr size_t dimension = 5;
   linear_system.dimension = dimension;
   linear_system.number_nonzeros = 7;
   // lower triangle, 0-based indexing
   linear_system.matrix_row_indices = {0, 1, 2, 4, 2, 3, 4};
   linear_system.matrix_column_indices = {0, 0, 1, 1, 2, 2, 4};
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

TEST(MA86Solver, Inertia) {
   MA86Solver solver;
   COOLinearSystem& linear_system = solver.get_coo_linear_system();
   linear_system.dimension = 5;
   linear_system.number_nonzeros = 7;
   linear_system.matrix_row_indices = {0, 1, 2, 4, 2, 3, 4};
   linear_system.matrix_column_indices = {0, 0, 1, 1, 2, 2, 4};
   linear_system.matrix_values = {2., 3., 4., 6., 1., 5., 1.};
   solver.initialize_memory();
   solver.do_symbolic_analysis();
   solver.do_numerical_factorization(false);
   const auto [number_positive, number_negative, number_zero] = solver.get_inertia();
   ASSERT_EQ(number_positive, 3);
   ASSERT_EQ(number_negative, 2);
   ASSERT_EQ(number_zero, 0);
}

TEST(MA86Solver, SingularMatrix) {
   // comes from hs015 solved with byrd preset (lower triangle, 0-based); relies on the summation of
   // duplicate (row, column) entries during the COO -> CSC conversion
   MA86Solver solver;
   COOLinearSystem& linear_system = solver.get_coo_linear_system();
   linear_system.dimension = 4;
   linear_system.number_nonzeros = 7;
   linear_system.matrix_row_indices = {0, 0, 1, 1, 1, 2, 3};
   linear_system.matrix_column_indices = {0, 0, 0, 1, 1, 2, 3};
   linear_system.matrix_values = {-0.0198, 0.625075, -0.277512, -0.624975, 0.625075, 0., 0.};
   solver.initialize_memory();
   solver.do_symbolic_analysis();
   solver.do_numerical_factorization(false);

   // expected inertia (1, 1, 2)
   ASSERT_TRUE(solver.matrix_is_singular());
}

TEST(MA86Solver, MultipleRHS) {
   MA86Solver solver;
   COOLinearSystem& linear_system = solver.get_coo_linear_system();
   constexpr size_t dimension = 5;
   linear_system.dimension = dimension;
   linear_system.number_nonzeros = 7;
   linear_system.rhs.resize(dimension);
   linear_system.matrix_row_indices = {0, 1, 2, 4, 2, 3, 4};
   linear_system.matrix_column_indices = {0, 0, 1, 1, 2, 2, 4};
   linear_system.matrix_values = {2., 3., 4., 6., 1., 5., 1.};
   solver.initialize_memory();
   solver.do_symbolic_analysis();
   solver.do_numerical_factorization(false);

   // two right-hand sides stored column-major: column 0 -> [1,2,3,4,5], column 1 -> [1,1,1,1,1]
   const std::array<double, 2 * dimension> rhs{8., 45., 31., 15., 17.,  5., 13., 10., 5., 7.};
   std::array<double, 2 * dimension> solution{};
   solver.solve_indefinite_system(rhs.data(), solution.data(), 2);

   const std::array<double, 2 * dimension> reference{1., 2., 3., 4., 5.,  1., 1., 1., 1., 1.};
   const double tolerance = 1e-8;
   for (size_t index: Range(2 * dimension)) {
      EXPECT_NEAR(solution[index], reference[index], tolerance);
   }
}

TEST(MA86Solver, SystemSize5FortranIndexing) {
   MA86Solver solver(Indexing::Fortran_indexing);
   COOLinearSystem& linear_system = solver.get_coo_linear_system();
   constexpr size_t dimension = 5;
   linear_system.dimension = dimension;
   linear_system.number_nonzeros = 7;
   // lower triangle, 1-based (Fortran) indexing
   linear_system.matrix_row_indices = {1, 2, 3, 5, 3, 4, 5};
   linear_system.matrix_column_indices = {1, 1, 2, 2, 3, 3, 5};
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

TEST(MA86Solver, InertiaFortranIndexing) {
   MA86Solver solver(Indexing::Fortran_indexing);
   COOLinearSystem& linear_system = solver.get_coo_linear_system();
   linear_system.dimension = 5;
   linear_system.number_nonzeros = 7;
   linear_system.matrix_row_indices = {1, 2, 3, 5, 3, 4, 5};
   linear_system.matrix_column_indices = {1, 1, 2, 2, 3, 3, 5};
   linear_system.matrix_values = {2., 3., 4., 6., 1., 5., 1.};
   solver.initialize_memory();
   solver.do_symbolic_analysis();
   solver.do_numerical_factorization(false);
   const auto [number_positive, number_negative, number_zero] = solver.get_inertia();
   ASSERT_EQ(number_positive, 3);
   ASSERT_EQ(number_negative, 2);
   ASSERT_EQ(number_zero, 0);
}

TEST(MA86Solver, SingularMatrixFortranIndexing) {
   // same singular system as SingularMatrix, in 1-based (Fortran) indexing
   MA86Solver solver(Indexing::Fortran_indexing);
   COOLinearSystem& linear_system = solver.get_coo_linear_system();
   linear_system.dimension = 4;
   linear_system.number_nonzeros = 7;
   linear_system.matrix_row_indices = {1, 1, 2, 2, 2, 3, 4};
   linear_system.matrix_column_indices = {1, 1, 1, 2, 2, 3, 4};
   linear_system.matrix_values = {-0.0198, 0.625075, -0.277512, -0.624975, 0.625075, 0., 0.};
   solver.initialize_memory();
   solver.do_symbolic_analysis();
   solver.do_numerical_factorization(false);

   // expected inertia (1, 1, 2)
   ASSERT_TRUE(solver.matrix_is_singular());
}

TEST(MA86Solver, MultipleRHSFortranIndexing) {
   MA86Solver solver(Indexing::Fortran_indexing);
   COOLinearSystem& linear_system = solver.get_coo_linear_system();
   constexpr size_t dimension = 5;
   linear_system.dimension = dimension;
   linear_system.number_nonzeros = 7;
   linear_system.rhs.resize(dimension);
   linear_system.matrix_row_indices = {1, 2, 3, 5, 3, 4, 5};
   linear_system.matrix_column_indices = {1, 1, 2, 2, 3, 3, 5};
   linear_system.matrix_values = {2., 3., 4., 6., 1., 5., 1.};
   solver.initialize_memory();
   solver.do_symbolic_analysis();
   solver.do_numerical_factorization(false);

   // two right-hand sides stored column-major: column 0 -> [1,2,3,4,5], column 1 -> [1,1,1,1,1]
   const std::array<double, 2 * dimension> rhs{8., 45., 31., 15., 17.,  5., 13., 10., 5., 7.};
   std::array<double, 2 * dimension> solution{};
   solver.solve_indefinite_system(rhs.data(), solution.data(), 2);

   const std::array<double, 2 * dimension> reference{1., 2., 3., 4., 5.,  1., 1., 1., 1., 1.};
   const double tolerance = 1e-8;
   for (size_t index: Range(2 * dimension)) {
      EXPECT_NEAR(solution[index], reference[index], tolerance);
   }
}
