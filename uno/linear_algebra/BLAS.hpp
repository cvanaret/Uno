// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BLAS_H
#define UNO_BLAS_H

#include "fortran_interface.h"
#define BLAS_copy_vector FC_GLOBAL_(dcopy, DCOPY)
#define BLAS_scale_vector FC_GLOBAL_(dscal, DSCAL)
#define BLAS_add_vectors FC_GLOBAL_(daxpy, DAXPY)
#define BLAS_dot_product FC_GLOBAL_(ddot, DDOT)
#define BLAS_triangular_back_solve FC_GLOBAL_(dtrsm, DTRSM)
#define BLAS_matrix_matrix_product FC_GLOBAL_(dgemm, DGEMM)
#define BLAS_triangular_matrix_matrix_product FC_GLOBAL_(dtrmm, DTRMM)

extern "C" {
   // y := x
   void BLAS_copy_vector(const int* n, const double* x, const int* incx, double* y, const int* incy);

   // x := alpha x
   void BLAS_scale_vector(const int* n, const double* alpha, double* x, const int* increment);

   // y := alpha x + y
   void BLAS_add_vectors(const int* n, const double* alpha, const double* x, const int* incx, double* y, const int* incy);

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

   // performs one of the matrix-matrix operations
   // C := alpha op(A) op(B) + beta C,
   // where
   // op(X) = X   or   op(X) = X^T
   void BLAS_matrix_matrix_product(const char* transa, const char* transb, const int* m, const int* n, const int* k,
      const double* alpha, const double* a, const int* lda, const double* b, const int* ldb, const double* beta,
      double* c, const int* ldc);

   void BLAS_triangular_matrix_matrix_product(char* side, char* uplo, char* transa, char* diag, int* m, int* n,
      double* alpha, double* a, int* lda, double* b, int* ldb);
}

#endif // UNO_BLAS_H