// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "linear_algebra/SymmetricMatrix.hpp"
#include "solvers/MA57/MA57Solver.hpp"

using namespace uno;

const size_t n = 5;
const size_t nnz = 7;
const std::array<double, n> reference{1., 2., 3., 4., 5.};

Vector<double> create_rhs() {
   Vector<double> rhs{8., 45., 31., 15., 17.};
   return rhs;
}

TEST(MA57Solver, SystemSize5) {
   SymmetricMatrix<size_t, double> matrix(n, nnz, false, "COO");
   matrix.insert(2., 0, 0);
   matrix.insert(3., 0, 1);
   matrix.insert(4., 1, 2);
   matrix.insert(6., 1, 4);
   matrix.insert(1., 2, 2);
   matrix.insert(5., 2, 3);
   matrix.insert(1., 4, 4);
   const Vector<double> rhs = create_rhs();
   Vector<double> result(n);
   result.fill(0.);

   MA57Solver solver(n, nnz);
   solver.do_symbolic_factorization(matrix);
   solver.do_numerical_factorization(matrix);
   solver.solve_indefinite_system(matrix, rhs, result);

   for (size_t index: Range(n)) {
      EXPECT_DOUBLE_EQ(result[index], reference[index]);
   }
}

TEST(MA57Solver, Inertia) {
   SymmetricMatrix<size_t, double> matrix(n, nnz, false, "COO");
   matrix.insert(2., 0, 0);
   matrix.insert(3., 0, 1);
   matrix.insert(4., 1, 2);
   matrix.insert(6., 1, 4);
   matrix.insert(1., 2, 2);
   matrix.insert(5., 2, 3);
   matrix.insert(1., 4, 4);
   const Vector<double> rhs = create_rhs();
   Vector<double> result(n);
   result.fill(0.);

   MA57Solver solver(n, nnz);
   solver.do_symbolic_factorization(matrix);
   solver.do_numerical_factorization(matrix);
   solver.solve_indefinite_system(matrix, rhs, result);

   const auto [number_positive, number_negative, number_zero] = solver.get_inertia();
   ASSERT_EQ(number_positive, 3);
   ASSERT_EQ(number_negative, 2);
   ASSERT_EQ(number_zero, 0);
}