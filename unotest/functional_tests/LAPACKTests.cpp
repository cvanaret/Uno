// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "fortran_interface.h"
#define LAPACK_cholesky_factorization FC_GLOBAL_(dpotrf, DPOTRF)

extern "C" {
   void LAPACK_cholesky_factorization(char* uplo, int* n, double* a, int* lda, int* info);
}

TEST(LAPACK, CholeskyFactorization) {
   int dimension = 3;
   // lower triangular matrix expressed in column-major order
   std::vector<double> matrix{
      4., 12., -16.,
      0., 37., -43.,
      0., 0., 98.
   };
   char uplo = 'L'; // lower triangular
   int info{0};
   LAPACK_cholesky_factorization(&uplo, &dimension, matrix.data(), &dimension, &info);


   ASSERT_EQ(info, 0);
   ASSERT_EQ(matrix[0], 2.);
   ASSERT_EQ(matrix[1], 6.);
   ASSERT_EQ(matrix[2], -8.);
   ASSERT_EQ(matrix[4], 1.);
   ASSERT_EQ(matrix[5], 5.);
   ASSERT_EQ(matrix[8], 3.);
}