// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LAPACK_H
#define UNO_LAPACK_H

#include "fortran_interface.h"
#define LAPACK_cholesky_factorization FC_GLOBAL_(dpotrf, DPOTRF)
#define LAPACK_symmetric_high_rank_update FC_GLOBAL_(dsyrk, DSYRK)

extern "C" {
   void LAPACK_cholesky_factorization(char* uplo, int* n, double* a, int* lda, int* info);

   // performs symmetric rank k update
   // C = alpha A A^T + beta C    or
   // C = alpha A^T A + beta C
   void LAPACK_symmetric_high_rank_update(char* uplo, char* trans, int* n, int* k, double* alpha, double* a, int* lda,
      double* beta, double* c, int* ldc);
}

#endif // UNO_LAPACK_H