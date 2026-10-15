// Copyright (c) 2018-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LINEARSOLVERTESTS_H
#define UNO_LINEARSOLVERTESTS_H

#include <array>
#include <cstddef>
#include <gtest/gtest.h>
#include "ingredients/subproblem_solvers/COOLinearSystem.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {

// number of times the pivoting tolerance can be increased from its default value (0: not supported).
// Specialized in the test file of each solver that supports it
template <typename SolverType>
constexpr size_t number_pivoting_tolerance_increases = 0;

// Type-parameterized test suite shared by every sparse indefinite linear solver interface (MA27, MA57, MA86, MUMPS, SSIDS)
template <typename SolverType>
class LinearSolverTest : public testing::Test {};

TYPED_TEST_SUITE_P(LinearSolverTest);

TYPED_TEST_P(LinearSolverTest, SystemSize5) {
   TypeParam solver;
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

TYPED_TEST_P(LinearSolverTest, Inertia) {
   TypeParam solver;
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

TYPED_TEST_P(LinearSolverTest, SingularMatrix) {
   // comes from hs015 solved with byrd preset
   TypeParam solver;
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

TYPED_TEST_P(LinearSolverTest, MultipleRHS) {
   TypeParam solver;
   COOLinearSystem& linear_system = solver.get_coo_linear_system();
   constexpr size_t dimension = 5;
   linear_system.dimension = dimension;
   linear_system.number_nonzeros = 7;
   linear_system.rhs.resize(dimension);
   // indices with Fortran-based indexing
   linear_system.matrix_row_indices = {1, 1, 2, 2, 3, 3, 5};
   linear_system.matrix_column_indices = {1, 2, 3, 5, 3, 4, 5};
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

TYPED_TEST_P(LinearSolverTest, MatrixVectorProduct) {
   TypeParam solver;
   COOLinearSystem& linear_system = solver.get_coo_linear_system();
   constexpr size_t dimension = 5;
   linear_system.dimension = dimension;
   linear_system.number_nonzeros = 8;
   linear_system.matrix_row_indices = {1, 1, 2, 2, 3, 3, 5, 3}; // duplicate (3, 3) entry
   linear_system.matrix_column_indices = {1, 2, 3, 5, 3, 4, 5, 3};
   linear_system.matrix_values = {2., 3., 4., 6., 0.5, 5., 1., 0.5};
   const std::array<double, dimension> vector{1., 2., 3., 4., 5.};
   std::array<double, dimension> result{};
   linear_system.compute_matrix_vector_product(vector.data(), result.data());

   const std::array<double, dimension> reference{8., 45., 31., 15., 17.};
   for (size_t index: Range(dimension)) {
      EXPECT_DOUBLE_EQ(result[index], reference[index]);
   }
}

TYPED_TEST_P(LinearSolverTest, PivotingToleranceIncrease) {
   TypeParam solver;
   COOLinearSystem& linear_system = solver.get_coo_linear_system();
   constexpr size_t dimension = 5;
   linear_system.dimension = dimension;
   linear_system.number_nonzeros = 7;
   linear_system.matrix_row_indices = {1, 1, 2, 2, 3, 3, 5};
   linear_system.matrix_column_indices = {1, 2, 3, 5, 3, 4, 5};
   linear_system.matrix_values = {2., 3., 4., 6., 1., 5., 1.};
   linear_system.rhs = {8., 45., 31., 15., 17.};
   solver.initialize_memory();
   solver.do_symbolic_analysis();
   solver.do_numerical_factorization(false);

   // the tolerance can be increased a finite number of times, after which the solver refuses
   size_t number_increases = 0;
   while (solver.increase_pivoting_tolerance()) {
      ++number_increases;
      ASSERT_LE(number_increases, number_pivoting_tolerance_increases<TypeParam>);
      // the refactorized matrix with the increased tolerance gives the same solution
      solver.do_numerical_factorization(false);
      Vector<double> result(dimension);
      solver.solve_indefinite_system(result.data());
      const std::array<double, dimension> reference{1., 2., 3., 4., 5.};
      for (size_t index: Range(dimension)) {
         EXPECT_NEAR(result[index], reference[index], 1e-8);
      }
   }
   EXPECT_EQ(number_increases, number_pivoting_tolerance_increases<TypeParam>);
   EXPECT_FALSE(solver.increase_pivoting_tolerance());

   // the increases can be undone, exactly as many times
   size_t number_decreases = 0;
   while (solver.decrease_pivoting_tolerance()) {
      ++number_decreases;
      ASSERT_LE(number_decreases, number_increases);
   }
   EXPECT_EQ(number_decreases, number_increases);
   // back at the default tolerance: the schedule can be replayed
   EXPECT_EQ(solver.increase_pivoting_tolerance(), 0 < number_increases);
}

REGISTER_TYPED_TEST_SUITE_P(LinearSolverTest,
   SystemSize5, Inertia, SingularMatrix, MultipleRHS, MatrixVectorProduct, PivotingToleranceIncrease);

} // namespace uno

#endif // UNO_LINEARSOLVERTESTS_H