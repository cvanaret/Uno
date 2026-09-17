// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include <array>
#include "ingredients/subproblem_solvers/HiPO/HiPOSolver.hpp"
#include "linear_algebra/Vector.hpp"

using namespace uno;

// HiPO is a linear solver for symmetric quasi-definite (SQD) systems, as they arise in interior-point
// methods: a positive-definite primal (1, 1) block and a negative-definite dual (2, 2) block. The
// solver only reads the lower triangle of the matrix (CSC, C-based indexing) and needs the expected
// inertia (set via set_expected_inertia) to set the expected sign of each pivot.
//
//   [[ 4, 1,  1 ],
//    [ 1, 3,  1 ],
//    [ 1, 1, -0.5]]    (2 primal variables, 1 dual constraint), inertia (2, 1, 0)

TEST(HiPOSolver, SystemSize3) {
   HiPOSolver solver;
   COOLinearSystem &linear_system = solver.get_coo_linear_system();
   constexpr size_t dimension = 3;
   linear_system.dimension = dimension;
   linear_system.number_nonzeros = 6;
   // lower triangle, C-based indexing
   linear_system.matrix_row_indices = {0, 1, 2, 1, 2, 2};
   linear_system.matrix_column_indices = {0, 0, 0, 1, 1, 2};
   linear_system.matrix_values = {4., 1., 1., 3., 1., -0.5};
   linear_system.rhs = {9., 10., 1.5};
   solver.initialize_memory();
   solver.set_expected_inertia(Inertia{2, 1, 0});
   solver.do_symbolic_analysis();
   solver.do_numerical_factorization(false);
   Vector<double> result(dimension);
   result.fill(0.);
   solver.solve_indefinite_system(result.data());

   const std::array<double, dimension> reference{1., 2., 3.};
   const double tolerance = 1e-8;
   for (size_t index: Range(dimension)) {
      EXPECT_NEAR(result[index], reference[index], tolerance);
   }
}

TEST(HiPOSolver, SPD) {
   // positive-definite (Hessian-only) regime: all pivots are expected positive (number_constraints = 0)
   //   [[4, 1, 1],
   //    [1, 3, 1],
   //    [1, 1, 3]]    inertia (3, 0, 0)
   HiPOSolver solver;
   COOLinearSystem &linear_system = solver.get_coo_linear_system();
   constexpr size_t dimension = 3;
   linear_system.dimension = dimension;
   linear_system.number_nonzeros = 6;
   linear_system.matrix_row_indices = {0, 1, 2, 1, 2, 2};
   linear_system.matrix_column_indices = {0, 0, 0, 1, 1, 2};
   linear_system.matrix_values = {4., 1., 1., 3., 1., 3.};
   linear_system.rhs = {6., 5., 5.};
   solver.initialize_memory();
   solver.set_expected_inertia(Inertia{3, 0, 0}); // SPD: all pivots expected positive
   solver.do_symbolic_analysis();
   solver.do_numerical_factorization(true);
   Vector<double> result(dimension);
   result.fill(0.);
   solver.solve_indefinite_system(result.data());

   const std::array<double, dimension> reference{1., 1., 1.};
   const double tolerance = 1e-8;
   for (size_t index: Range(dimension)) {
      EXPECT_NEAR(result[index], reference[index], tolerance);
   }
   const auto [number_positive, number_negative, number_zero] = solver.get_inertia();
   ASSERT_EQ(number_positive, 3);
   ASSERT_EQ(number_negative, 0);
   ASSERT_EQ(number_zero, 0);
}

TEST(HiPOSolver, Inertia) {
   HiPOSolver solver;
   COOLinearSystem &linear_system = solver.get_coo_linear_system();
   linear_system.dimension = 3;
   linear_system.number_nonzeros = 6;
   linear_system.matrix_row_indices = {0, 1, 2, 1, 2, 2};
   linear_system.matrix_column_indices = {0, 0, 0, 1, 1, 2};
   linear_system.matrix_values = {4., 1., 1., 3., 1., -0.5};
   solver.initialize_memory();
   solver.set_expected_inertia(Inertia{2, 1, 0});
   solver.do_symbolic_analysis();
   solver.do_numerical_factorization(false);
   const auto [number_positive, number_negative, number_zero] = solver.get_inertia();
   ASSERT_EQ(number_positive, 2);
   ASSERT_EQ(number_negative, 1);
   ASSERT_EQ(number_zero, 0);
   ASSERT_FALSE(solver.matrix_is_singular());
}

TEST(HiPOSolver, DuplicateEntries) {
   // the same matrix, but with the (0, 0) and (2, 2) entries each split into two duplicate COO entries
   // that the solver must sum: 4 = 1.5 + 2.5 and -0.5 = -0.2 + (-0.3)
   HiPOSolver solver;
   COOLinearSystem &linear_system = solver.get_coo_linear_system();
   constexpr size_t dimension = 3;
   linear_system.dimension = dimension;
   linear_system.number_nonzeros = 8;
   linear_system.matrix_row_indices = {0, 0, 1, 2, 1, 2, 2, 2};
   linear_system.matrix_column_indices = {0, 0, 0, 0, 1, 1, 2, 2};
   linear_system.matrix_values = {1.5, 2.5, 1., 1., 3., 1., -0.2, -0.3};
   linear_system.rhs = {9., 10., 1.5};
   solver.initialize_memory();
   solver.set_expected_inertia(Inertia{2, 1, 0});
   solver.do_symbolic_analysis();
   solver.do_numerical_factorization(false);
   Vector<double> result(dimension);
   result.fill(0.);
   solver.solve_indefinite_system(result.data());

   const std::array<double, dimension> reference{1., 2., 3.};
   const double tolerance = 1e-8;
   for (size_t index: Range(dimension)) {
      EXPECT_NEAR(result[index], reference[index], tolerance);
   }
}

TEST(HiPOSolver, MultipleRHS) {
   // native multiple right-hand side solve on the SQD matrix: two columns solved in one call.
   //   column 0: rhs [9, 10, 1.5] -> x [1, 2, 3]
   //   column 1: rhs [9,  3, 1.5] -> x [2, 0, 1]
   HiPOSolver solver;
   COOLinearSystem &linear_system = solver.get_coo_linear_system();
   linear_system.dimension = 3;
   linear_system.number_nonzeros = 6;
   linear_system.matrix_row_indices = {0, 1, 2, 1, 2, 2};
   linear_system.matrix_column_indices = {0, 0, 0, 1, 1, 2};
   linear_system.matrix_values = {4., 1., 1., 3., 1., -0.5};
   solver.initialize_memory();
   solver.set_expected_inertia(Inertia{2, 1, 0});
   solver.do_symbolic_analysis();
   solver.do_numerical_factorization(false);

   // right-hand sides and solutions are column-major blocks (2 columns of length 3)
   const std::array<double, 6> rhs{9., 10., 1.5, 9., 3., 1.5};
   std::array<double, 6> solution{};
   solver.solve_indefinite_system(rhs.data(), solution.data(), 2);

   const std::array<double, 6> reference{1., 2., 3., 2., 0., 1.};
   const double tolerance = 1e-8;
   for (size_t index: Range(6)) {
      EXPECT_NEAR(solution[index], reference[index], tolerance);
   }
}

// General symmetric indefinite saddle-point matrix (the one used by the MA57/SSIDS tests), nonsingular
// with inertia (3, 2, 0). It has missing diagonal entries (e.g. node 3 has no diagonal, so in CSC its
// lower-triangle column would be empty) - the shape of an augmented system with a zero (2, 2) block and
// no dual regularization. HiPO requires an explicit diagonal in every column, so HiPOSolver inserts the
// missing diagonals (with value 0) during the COO -> CSC conversion; the factorization is then correct.
TEST(HiPOSolver, SymmetricIndefinite5x5) {
   HiPOSolver solver;
   COOLinearSystem &linear_system = solver.get_coo_linear_system();
   constexpr size_t dimension = 5;
   linear_system.dimension = dimension;
   linear_system.number_nonzeros = 7;
   // lower triangle, C-based indexing
   linear_system.matrix_row_indices = {0, 1, 2, 4, 2, 3, 4};
   linear_system.matrix_column_indices = {0, 0, 1, 1, 2, 2, 4};
   linear_system.matrix_values = {2., 3., 4., 6., 1., 5., 1.};
   linear_system.rhs = {8., 45., 31., 15., 17.};
   solver.initialize_memory();
   solver.set_expected_inertia(Inertia{3, 2, 0});
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
   const auto [number_positive, number_negative, number_zero] = solver.get_inertia();
   ASSERT_EQ(number_positive, 3);
   ASSERT_EQ(number_negative, 2);
   ASSERT_EQ(number_zero, 0);
}

TEST(HiPOSolver, SingularMatrix) {
   // rank-deficient 4x4 matrix: the trailing block [[1, 1], [1, 1]] is singular, so the matrix has
   // inertia (3, 0, 1) (one zero eigenvalue). Checks that HiPO reports the zero pivot
   // (get_inertia / matrix_is_singular) rather than silently regularizing it.
   //   [[2, 1, 0, 0],
   //    [1, 2, 0, 0],
   //    [0, 0, 1, 1],
   //    [0, 0, 1, 1]]
   HiPOSolver solver;
   COOLinearSystem &linear_system = solver.get_coo_linear_system();
   constexpr size_t dimension = 4;
   linear_system.dimension = dimension;
   linear_system.number_nonzeros = 6;
   linear_system.matrix_row_indices = {0, 1, 1, 2, 3, 3};
   linear_system.matrix_column_indices = {0, 0, 1, 2, 2, 3};
   linear_system.matrix_values = {2., 1., 2., 1., 1., 1.};
   solver.initialize_memory();
   solver.set_expected_inertia(Inertia{3, 0, 1});
   solver.do_symbolic_analysis();
   solver.do_numerical_factorization(false);

   ASSERT_TRUE(solver.matrix_is_singular());
   const auto [number_positive, number_negative, number_zero] = solver.get_inertia();
   ASSERT_EQ(number_positive, 3);
   ASSERT_EQ(number_negative, 0);
   ASSERT_EQ(number_zero, 1);
}
