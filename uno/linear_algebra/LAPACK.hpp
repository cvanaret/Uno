// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LAPACK_H
#define UNO_LAPACK_H

#include "fortran_interface.h"
#define LAPACK_cholesky_factorization FC_GLOBAL_(dpotrf, DPOTRF)
#define LAPACK_symmetric_high_rank_update FC_GLOBAL_(dsyrk, DSYRK)
#define LAPACK_triangular_back_solve FC_GLOBAL_(dtrsm, DTRSM)
//#define LAPACK_matrix_matrix_product FC_GLOBAL_(dgemm, DGEMM)
#define LAPACK_triangular_matrix_matrix_product FC_GLOBAL_(dtrmm, DTRMM)
#define LAPACK_add_vector FC_GLOBAL_(daxpy, DAXPY)

extern "C" {
   // perform Cholesky factorization of A
   // A = U^T U    or
   // A = L L^T
   void LAPACK_cholesky_factorization(char* uplo, int* n, double* a, int* lda, int* info);

   // performs symmetric rank k update:
   // C = alpha A A^T + beta C    or
   // C = alpha A^T A + beta C
   void LAPACK_symmetric_high_rank_update(char* uplo, char* trans, int* n, int* k, double* alpha, double* a, int* lda,
      double* beta, double* c, int* ldc);

   // performs back-solve with lower-triangular A:
   // op(A) X = alpha B    or
   // X op(A) = alpha B
   // where
   // op(A) = A   or   op(A) = A^T
   // X is overwritten on B
   void LAPACK_triangular_back_solve(char* side, char* uplo, char* transa, char* diag, int* m, int* n, double* alpha,
      double* a, int* lda, double* b, int* ldb);

   /*
   void LAPACK_matrix_matrix_product(char* transa, char* transb, int* m, int* n, int* k, double* alpha, double* a,
      int* lda, double* b, int* ldb, double* beta, double* c, int* ldc);
   */

   void LAPACK_triangular_matrix_matrix_product(char* side, char* uplo, char* transa, char* diag, int* m, int* n,
      double* alpha, double* a, int* lda, double* b, int* ldb);

   // y := alpha x + y
   void LAPACK_add_vector(const int* n, const double* alpha, const double* x, const int* incx, const double* y, int* incy);
}

#endif // UNO_LAPACK_H