// Copyright (c) 2025-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include <cmath>
#include "linear_algebra/BLAS.hpp"

TEST(BLAS, DotProduct) {
   constexpr size_t dimension = 3;
   // lower triangular matrix expressed in column-major order
   std::vector<double> x{-1., 2., -3.};
   std::vector<double> y{4., -5., 6.};
   const double dot_product = uno::blas1::dot(dimension, x.data(), y.data());
   ASSERT_EQ(dot_product, -32.);
}

TEST(BLAS, TriangularBackSolve) {
   constexpr size_t dimension = 3;
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
   uno::blas3::triangular_back_solve('L', 'L', 'N', 'N', dimension, dimension, 1., matrix.data(), dimension, rhs.data(), dimension);
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

/*
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
   uno::blas3::tri
   BLAS_triangular_matrix_matrix_product(&side, &uplo, &transa, &diag, &m, &n, &alpha, lower_triangular_matrix.data(),
      &lda, matrix.data(), &ldb);
   ASSERT_EQ(matrix[0], 8.);
   ASSERT_EQ(matrix[1], -13.);
   ASSERT_EQ(matrix[2], 795.);
   ASSERT_EQ(matrix[3], 12.);
   ASSERT_EQ(matrix[4], 295.);
   ASSERT_EQ(matrix[5], -545.);
}
*/

// compute high-rank update I + v v^T
TEST(BLAS, SymmetricHighRank1Update) {
   constexpr size_t dimension = 3;
   // identity matrix
   std::vector<double> matrix{
      1., 0., 0.,
      0., 1., 0.,
      0., 0., 1.
   };
   std::vector<double> vector{2., 3., 4.};
   uno::blas3::symmetric_rank_k_update('L', 'N', dimension, dimension, 1, 1., vector.data(), dimension, 1.,
      matrix.data(), dimension);
   // output matrix should be
   //   5., 6., 8.,
   //   0., 10., 12.,
   //   0., 0., 17.
   ASSERT_EQ(matrix[0], 5.);
   ASSERT_EQ(matrix[1], 6.);
   ASSERT_EQ(matrix[2], 8.);
   ASSERT_EQ(matrix[4], 10.);
   ASSERT_EQ(matrix[5], 12.);
   ASSERT_EQ(matrix[8], 17.);
}

// compute high-rank update 0 + L L^T
TEST(BLAS, SymmetricHighRank3Update) {
   constexpr size_t dimension = 3;
   std::vector<double> matrix{
      0., 0., 0.,
      0., 0., 0.,
      0., 0., 0.
   };
   const double invsqrt10 = 1./std::sqrt(10.);
   std::vector<double> lower_triangular_matrix{
      0., invsqrt10, 2.*invsqrt10,
      0., 0., 3.*invsqrt10,
      0., 0., 0.
   };
   uno::blas3::symmetric_rank_k_update('L', 'N', dimension, dimension, dimension, 1., lower_triangular_matrix.data(),
      dimension, 1., matrix.data(), dimension);
   constexpr double tolerance = 1e-6;
   // output matrix should be
   //    0., 0., 0.,
   //    0., 0.1, 0.2,
   //    0., 0., 1.3
   EXPECT_NEAR(matrix[4], 0.1, tolerance);
   EXPECT_NEAR(matrix[5], 0.2, tolerance);
   EXPECT_NEAR(matrix[8], 1.3, tolerance);
}