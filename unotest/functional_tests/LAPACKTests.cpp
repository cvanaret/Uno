// Copyright (c) 2025-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "linear_algebra/LAPACK.hpp"
#include "linear_algebra/LAPACK_extension.hpp"

TEST(LAPACK, CholeskyFactorization) {
   constexpr size_t dimension = 3;
   // lower triangular matrix expressed in column-major order
   std::vector<double> matrix{
      4., 12., -16.,
      0., 37., -43.,
      0., 0., 98.
   };
   const bool success = uno::lapack::cholesky_factorization('L', dimension, matrix.data(), dimension);

   ASSERT_EQ(success, true);
   ASSERT_EQ(matrix[0], 2.);
   ASSERT_EQ(matrix[1], 6.);
   ASSERT_EQ(matrix[2], -8.);
   ASSERT_EQ(matrix[4], 1.);
   ASSERT_EQ(matrix[5], 5.);
   ASSERT_EQ(matrix[8], 3.);
}

TEST(LAPACK, LDLtFactorization2) {
   constexpr size_t dimension = 2;
   // lower triangular matrix expressed in column-major order
   std::vector<double> matrix{
      3., -6.,
      0., 1.,
   };
   const bool success = uno::ldlt_nopiv_lvl2_rightlooking(matrix.data(), dimension, dimension, 0.);

   ASSERT_EQ(success, true);
   ASSERT_EQ(matrix[0], 3.);
   ASSERT_EQ(matrix[1], -2.);
   ASSERT_EQ(matrix[3], -11.);
}

TEST(LAPACK, LDLtFactorization2Singular) {
   constexpr size_t dimension = 2;
   // lower triangular matrix expressed in column-major order
   std::vector<double> matrix{
      3., -6.,
      0., 12.,
   };
   const bool success = uno::ldlt_nopiv_lvl2_rightlooking(matrix.data(), dimension, dimension, 0.);
   // the matrix is singular, expect a failure
   ASSERT_EQ(success, false);
}

TEST(LAPACK, LDLtFactorization3) {
   constexpr size_t dimension = 3;
   // lower triangular matrix expressed in column-major order
   std::vector<double> matrix{
      4., 12., -16.,
      0., 37., -43.,
      0., 0., 98.
   };
   const bool success = uno::ldlt_nopiv_lvl2_rightlooking(matrix.data(), dimension, dimension, 0.);

   ASSERT_EQ(success, true);
   ASSERT_EQ(matrix[0], 4.);
   ASSERT_EQ(matrix[1], 3.);
   ASSERT_EQ(matrix[2], -4.);
   ASSERT_EQ(matrix[4], 1.);
   ASSERT_EQ(matrix[5], 5.);
   ASSERT_EQ(matrix[8], 9.);
}