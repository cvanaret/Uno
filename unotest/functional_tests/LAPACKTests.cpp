// Copyright (c) 2025-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "linear_algebra/LAPACK.hpp"

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