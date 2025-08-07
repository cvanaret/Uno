// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "ingredients/subproblem_solvers/MA27/MA27Solver.hpp"

using namespace uno;

/*
TEST(MA27Solver, SystemSize5) {
   const size_t n = 5;
   const size_t nnz = 7;
   SymmetricMatrix<size_t, double> matrix("COO", n, nnz, 0);
   matrix.insert(2., 0, 0);
   matrix.insert(3., 0, 1);
   matrix.insert(4., 1, 2);
   matrix.insert(6., 1, 4);
   matrix.insert(1., 2, 2);
   matrix.insert(5., 2, 3);
   matrix.insert(1., 4, 4);
   const Vector<double> rhs{8., 45., 31., 15., 17.};
   Vector<double> result(n);
   result.fill(0.);
   const std::array<double, n> reference{1., 2., 3., 4., 5.};

   MA27Solver solver;
   solver.initialize_memory(n, 0, nnz, 0);
   solver.do_symbolic_analysis(matrix);
   solver.do_numerical_factorization(matrix);
   solver.solve_indefinite_system(matrix, rhs, result);

   for (size_t index: Range(n)) {
      EXPECT_DOUBLE_EQ(result[index], reference[index]);
   }
}

TEST(MA27Solver, Inertia) {
   const size_t n = 5;
   const size_t nnz = 7;
   SymmetricMatrix<size_t, double> matrix("COO", n, nnz, 0);
   matrix.insert(2., 0, 0);
   matrix.insert(3., 0, 1);
   matrix.insert(4., 1, 2);
   matrix.insert(6., 1, 4);
   matrix.insert(1., 2, 2);
   matrix.insert(5., 2, 3);
   matrix.insert(1., 4, 4);

   MA27Solver solver;
   solver.initialize_memory(n, 0, nnz, 0);
   solver.do_symbolic_analysis(matrix);
   solver.do_numerical_factorization(matrix);

   const auto [number_positive, number_negative, number_zero] = solver.get_inertia();
   ASSERT_EQ(number_positive, 3);
   ASSERT_EQ(number_negative, 2);
   ASSERT_EQ(number_zero, 0);
}

TEST(MA27Solver, SingularMatrix) {
   const size_t n = 4;
   const size_t nnz = 7;
   // comes from hs015 solved with byrd preset
   SymmetricMatrix<size_t, double> matrix("COO", n, nnz, 0);
   matrix.insert(-0.0198, 0, 0);
   matrix.insert(0.625075, 0, 0);
   matrix.insert(-0.277512, 0, 1);
   matrix.insert(-0.624975, 1, 1);
   matrix.insert(0.625075, 1, 1);
   matrix.insert(0., 2, 2);
   matrix.insert(0., 3, 3);
   MA27Solver solver;
   solver.initialize_memory(n, 0, nnz, 0);
   solver.do_symbolic_analysis(matrix);
   solver.do_numerical_factorization(matrix);

   // expected inertia (1, 1, 2)
   ASSERT_TRUE(solver.matrix_is_singular());
}
*/