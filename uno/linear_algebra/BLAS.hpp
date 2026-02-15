// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BLAS_H
#define UNO_BLAS_H

#include "fortran_interface.h"
#define BLAS_add_vector FC_GLOBAL_(daxpy, DAXPY)
#define BLAS_dot_product FC_GLOBAL_(ddot, DDOT)
#define BLAS_triangular_back_solve FC_GLOBAL_(dtrsm, DTRSM)
//#define BLAS_matrix_matrix_product FC_GLOBAL_(dgemm, DGEMM)
#define BLAS_triangular_matrix_matrix_product FC_GLOBAL_(dtrmm, DTRMM)

extern "C" {
   // y := alpha x + y
   void BLAS_add_vector(const int* n, const double* alpha, const double* x, const int* incx, const double* y, int* incy);

   // x^T y
   double BLAS_dot_product(const int* n, const double* x, const int* incx, const double* y, const int* incy);

   // performs back-solve with lower-triangular A:
   // op(A) X = alpha B    or
   // X op(A) = alpha B
   // where
   // op(A) = A   or   op(A) = A^T
   // X is overwritten on B
   void BLAS_triangular_back_solve(char* side, char* uplo, char* transa, char* diag, int* m, int* n, double* alpha,
      double* a, int* lda, double* b, int* ldb);

   /*
   void BLAS_matrix_matrix_product(char* transa, char* transb, int* m, int* n, int* k, double* alpha, double* a,
      int* lda, double* b, int* ldb, double* beta, double* c, int* ldc);
   */

   void BLAS_triangular_matrix_matrix_product(char* side, char* uplo, char* transa, char* diag, int* m, int* n,
      double* alpha, double* a, int* lda, double* b, int* ldb);
}

#endif // UNO_BLAS_H