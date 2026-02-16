// Copyright (c) 2025-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "linear_algebra/LAPACK.hpp"

TEST(LAPACK, CholeskyFactorization) {
   int dimension = 3;
   // lower triangular matrix expressed in column-major order
   std::vector<double> matrix{
      4., 12., -16.,
      0., 37., -43.,
      0., 0., 98.
   };
   char uplo = 'L'; // lower triangular
   int info = 0;
   LAPACK_cholesky_factorization(&uplo, &dimension, matrix.data(), &dimension, &info);

   ASSERT_EQ(info, 0);
   ASSERT_EQ(matrix[0], 2.);
   ASSERT_EQ(matrix[1], 6.);
   ASSERT_EQ(matrix[2], -8.);
   ASSERT_EQ(matrix[4], 1.);
   ASSERT_EQ(matrix[5], 5.);
   ASSERT_EQ(matrix[8], 3.);
}

// compute high-rank update I + v v^T
TEST(LAPACK, SymmetricHighRank1Update) {
   size_t dimension = 3;
   // identity matrix
   std::vector<double> matrix{
      1., 0., 0.,
      0., 1., 0.,
      0., 0., 1.
   };
   std::vector<double> vector{2., 3., 4.};
   char uplo = 'L'; // lower triangular
   char trans = 'N';
   int dim = static_cast<int>(dimension);
   int k = 1; // rank of the update
   double alpha = 1.;
   double beta = 1.;
   LAPACK_symmetric_high_rank_update(&uplo, &trans, &dim, &k, &alpha, vector.data(), &dim, &beta,
      matrix.data(), &dim);
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
TEST(LAPACK, SymmetricHighRank3Update) {
   size_t dimension = 3;
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
   char uplo = 'L'; // lower triangular
   char trans = 'N';
   int dim = static_cast<int>(dimension);
   int k = dim; // rank of the update
   double alpha = 1.;
   double beta = 1.;
   LAPACK_symmetric_high_rank_update(&uplo, &trans, &dim, &k, &alpha, lower_triangular_matrix.data(), &dim, &beta,
      matrix.data(), &dim);
   const double tolerance = 1e-6;
   // output matrix should be
   //    0., 0., 0.,
   //    0., 0.1, 0.2,
   //    0., 0., 1.3
   EXPECT_NEAR(matrix[4], 0.1, tolerance);
   EXPECT_NEAR(matrix[5], 0.2, tolerance);
   EXPECT_NEAR(matrix[8], 1.3, tolerance);
}