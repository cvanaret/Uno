// Copyright (c) 2025-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include <cmath>
#include "linear_algebra/BLAS.hpp"

TEST(BLAS, DotProduct) {
   const int dimension = 3;
   // lower triangular matrix expressed in column-major order
   std::vector<double> x{-1., 2., -3.};
   std::vector<double> y{4., -5., 6.};
   const int increment = 1;
   const double dot_produt = BLAS_dot_product(&dimension, x.data(), &increment, y.data(), &increment);
   ASSERT_EQ(dot_produt, -32.);
}

TEST(BLAS, TriangularBackSolve) {
   int dimension = 3;
   // lower triangular matrix expressed in column-major order
   std::vector<double> matrix{
      1., -2., 3.5,
      0., 1., -2.5,
      0., 0., 1.
   };
   std::vector<double> rhs{
      1., 0., 0,
      0., 1., 0,
      0., 0., 1.
   };
   // solve A X = alpha B
   char side = 'L'; // left
   char uplo = 'L'; // lower triangular
   char transa = 'N'; // not transposed
   char diag = 'N'; // not unit triangular
   int m = dimension; // number of rows of B
   int n = dimension; // number of RHS (columns of B)
   double alpha = 1.;
   int lda = dimension; // leading dimension of A
   int ldb = dimension; // leading dimension of B
   BLAS_triangular_back_solve(&side, &uplo, &transa, &diag, &m, &n, &alpha, matrix.data(), &lda, rhs.data(), &ldb);
   ASSERT_EQ(rhs[0], 1.);
   ASSERT_EQ(rhs[1], 2.);
   ASSERT_EQ(rhs[2], 1.5);
   ASSERT_EQ(rhs[3], 0.);
   ASSERT_EQ(rhs[4], 1.);
   ASSERT_EQ(rhs[5], 2.5);
   ASSERT_EQ(rhs[6], 0.);
   ASSERT_EQ(rhs[7], 0.);
   ASSERT_EQ(rhs[8], 1.);
}

// B := alpha A B with lower triangular A
TEST(BLAS, TriangularMatrixMatrixProduct) {
   // lower triangular matrix expressed in column-major order
   std::vector<double> lower_triangular_matrix{ // 3x3
      4., 12., -16.,
      0., 37., -43.,
      0., 0., 98.
   };
   std::vector<double> matrix{ // 3x2
      2., -1., 8.,
      3., 7., -2.
   };

   char side = 'L';
   char uplo = 'L';
   char transa = 'N';
   char diag = 'N';
   int m = 3; // number of rows of matrix
   int n = 2; // number of columns of matrix
   double alpha = 1.;
   int lda = 3; // leading dimension of lower_triangular_matrix
   int ldb = 3; // leading dimension of matrix
   BLAS_triangular_matrix_matrix_product(&side, &uplo, &transa, &diag, &m, &n, &alpha, lower_triangular_matrix.data(),
      &lda, matrix.data(), &ldb);
   ASSERT_EQ(matrix[0], 8.);
   ASSERT_EQ(matrix[1], -13.);
   ASSERT_EQ(matrix[2], 795.);
   ASSERT_EQ(matrix[3], 12.);
   ASSERT_EQ(matrix[4], 295.);
   ASSERT_EQ(matrix[5], -545.);
}