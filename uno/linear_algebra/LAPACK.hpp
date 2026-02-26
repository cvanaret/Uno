// Copyright (c) 2025-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LAPACK_H
#define UNO_LAPACK_H

#include "fortran_interface.h"
#define LAPACK_cholesky_factorization FC_GLOBAL_(dpotrf, DPOTRF)
#define LAPACK_symmetric_high_rank_update FC_GLOBAL_(dsyrk, DSYRK)

extern "C" {
   // perform Cholesky factorization of A
   // A = U^T U    or
   // A = L L^T
   void LAPACK_cholesky_factorization(const char* uplo, const int* n, double* a, const int* lda, int* info);

   // performs symmetric rank k update:
   // C = alpha A A^T + beta C    or
   // C = alpha A^T A + beta C
   void LAPACK_symmetric_high_rank_update(const char* uplo, const char* trans, const int* n, const int* k, const double* alpha,
      const double* a, const int* lda, const double* beta, double* c, const int* ldc);
}

#endif // UNO_LAPACK_H