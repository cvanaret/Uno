// Copyright (c) 2025-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_LAPACK_H
#define UNO_LAPACK_H

#include <cassert>
#include "fortran_interface.h"
#define dpotrf FC_GLOBAL_(dpotrf, DPOTRF)
#define dsytrf FC_GLOBAL_(dsytrf, DSYTRF)
#define dsytrs FC_GLOBAL_(dsytrs, DSYTRS)

extern "C" {
   // performs Cholesky factorization of a symmetric positive definite matrix A
   // A = U^T U    or
   // A = L L^T
   void dpotrf(const char* uplo, const int* n, double* a, const int* lda, int* info);

   // performs the factorization of a symmetric matrix A using the Bunch-Kaufman diagonal pivoting method
   // A = U^T D U  or
   // A = L D L^T
   void dsytrf(const char* uplo, const int* n, double* a, const int* lda, int* ipiv, double* work, const int* lwork, int* info);

   // solves a system of linear equations A X = B with a symmetric matrix A using the factorization computed by dsytrf
   void dsytrs(const char* uplo, const int* n, const int* nrhs, const double* a, const int* lda, const int* ipiv, double* b,
      const int* ldb, int* info);
}

namespace uno {
   namespace lapack {
      // performs Cholesky factorization of a symmetric positive definite matrix A
      // A = U^T U    or
      // A = L L^T
      // returns true upon success, false upon failure
      inline bool cholesky_factorization(char uplo, size_t dimension, double* a, size_t leading_dimension) {
         const int n = static_cast<int>(dimension);
         const int lda = static_cast<int>(leading_dimension);
         assert(lda >= std::max(1, n));
         int info = 0;
         dpotrf(&uplo, &n, a, &lda, &info);
         return (info == 0);
      }

      // performs the factorization of a symmetric matrix A using the Bunch-Kaufman diagonal pivoting method
      // A = U^T D U  or
      // A = L D L^T
      inline std::pair<bool, std::vector<int>> bunch_kaufman_factorization(char uplo, size_t dimension, double* a,
            size_t leading_dimension) {
         const int n = static_cast<int>(dimension);
         const int lda = static_cast<int>(leading_dimension);
         assert(lda >= std::max(1, n));
         std::vector<int> ipiv(dimension);
         // first call to get the optimal lwork
         double work_size = 0.;
         int lwork = -1;
         int info = 0;
         dsytrf(&uplo, &n, a, &lda, ipiv.data(), &work_size, &lwork, &info);
         // second call to factorize
         lwork = static_cast<int>(work_size);
         assert(lwork >= 0);
         std::vector<double> work(static_cast<size_t>(lwork));
         dsytrf(&uplo, &n, a, &lda, ipiv.data(), work.data(), &lwork, &info);
         return {(info == 0), std::move(ipiv)};
      }

      inline bool bunch_kaufman_solve(char uplo, size_t dimension, const double* a, size_t leading_dimension,
            const int* ipiv, double* b) {
         const int n = static_cast<int>(dimension);
         constexpr int nrhs = 1;
         const int lda = static_cast<int>(leading_dimension);
         const int ldb = static_cast<int>(dimension);
         int info = 0;
         dsytrs(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
         return (info == 0);
      }
   }
} // namespace

#endif // UNO_LAPACK_H