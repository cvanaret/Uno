// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BLAS_H
#define UNO_BLAS_H

#include <cassert>
#include "fortran_interface.h"
#define dcopy FC_GLOBAL_(dcopy, DCOPY)
#define dscal FC_GLOBAL_(dscal, DSCAL)
#define daxpy FC_GLOBAL_(daxpy, DAXPY)
#define ddot FC_GLOBAL_(ddot, DDOT)
#define dgemv FC_GLOBAL_(dgemv, DGEMV)
#define dtrsm FC_GLOBAL_(dtrsm, DTRSM)
#define dgemm FC_GLOBAL_(dgemm, DGEMM)
#define dtrmm FC_GLOBAL_(dtrmm, DTRMM)
#define dsyr FC_GLOBAL_(dsyr, DSYR)
#define dsyrk FC_GLOBAL_(dsyrk, DSYRK)

extern "C" {
   // y := x
   void dcopy(const int* n, const double* x, const int* incx, double* y, const int* incy);

   // x := alpha x
   void dscal(const int* n, const double* alpha, double* x, const int* increment);

   // y := alpha x + y
   void daxpy(const int* n, const double* alpha, const double* x, const int* incx, double* y, const int* incy);

   // x^T y
   double ddot(const int* n, const double* x, const int* incx, const double* y, const int* incy);

   // performs the matrix-vector operations
   // y := alpha A x + beta y,   or   y := alpha A^T x + beta y,
   void dgemv(const char* trans, const int* m, const int* n, const double* alpha, const double* a,
      const int* lda, const double* x, const int* incx, const double* beta, double* y, const int* incy);

   // performs back-solve with lower-triangular A:
   // op(A) X = alpha B    or
   // X op(A) = alpha B
   // where
   // op(A) = A   or   op(A) = A^T
   // X is overwritten on B
   void dtrsm(const char* side, const char* uplo, const char* transa, const char* diag, const int* m,
      const int* n, const double* alpha, const double* a, const int* lda, double* b, const int* ldb);

   // performs one of the matrix-matrix operations
   // C := alpha op(A) op(B) + beta C,
   // where
   // op(X) = X   or   op(X) = X^T
   void dgemm(const char* transa, const char* transb, const int* m, const int* n, const int* k, const double* alpha,
      const double* a, const int* lda, const double* b, const int* ldb, const double* beta, double* c, const int* ldc);

   // performs one of the matrix-matrix operations
   // B := alpha op(A) B   or   B := alpha B op(A)
   // where
   // alpha is a scalar, B is an m by n matrix, A is a unit, or non-unit, triangular matrix, and
   // op(X) = X   or   op(X) = X^T
   void dtrmm(const char* side, const char* uplo, const char* transa, const char* diag, const int* m, const int* n,
      const double* alpha, const double* a, const int* lda, double* b, const int* ldb);

   // performs symmetric rank-1 update:
   // A := alpha x x^T + A
   void dsyr(const char* uplo, const int* n, const double* alpha, const double* x, const int* incx, double* a, const int* lda);

   // performs symmetric rank k update:
   // C = alpha A A^T + beta C    or
   // C = alpha A^T A + beta C
   void dsyrk(const char* uplo, const char* trans, const int* n, const int* k, const double* alpha, const double* a,
      const int* lda, const double* beta, double* c, const int* ldc);
}

namespace uno {
   namespace blas1 {
      // y := x
      inline void copy(size_t size, const double* x, double* y) {
         const int n = static_cast<int>(size);
         constexpr int increment = 1;
         dcopy(&n, x, &increment, y, &increment);
      }

      // x := alpha x
      inline void scale(size_t size, double alpha, double* x) {
         const int n = static_cast<int>(size);
         constexpr int increment = 1;
         dscal(&n, &alpha, x, &increment);
      }

      // y := alpha x + y
      inline void add(size_t size, double alpha, const double* x, double* y) {
         const int n = static_cast<int>(size);
         constexpr int increment = 1;
         daxpy(&n, &alpha, x, &increment, y, &increment);
      }

      // x^T y
      inline double dot(size_t size, const double* x, const double* y) {
         const int n = static_cast<int>(size);
         constexpr int increment = 1;
         return ddot(&n, x, &increment, y, &increment);
      }
   }

   namespace blas2 {
      // performs the matrix-vector operations
      // y := alpha A x + beta y,   or   y := alpha A^T x + beta y,
      inline void matrix_vector_product(char trans, size_t number_rows, size_t number_colums, double alpha, const double* a,
            size_t leading_dimension, const double* x, double beta, double* y) {
         const int m = static_cast<int>(number_rows);
         const int n = static_cast<int>(number_colums);
         const int lda = static_cast<int>(leading_dimension);
         assert(lda >= std::max(1, m));
         constexpr int increment = 1;
         dgemv(&trans, &m, &n, &alpha, a, &lda, x, &increment, &beta, y, &increment);
      }
   }

   namespace blas3 {
      // performs one of the matrix-matrix operations
      // C := alpha op(A) op(B) + beta C,
      // where
      // op(X) = X   or   op(X) = X^T
      inline void matrix_matrix_product(char transa, char transb, size_t number_rows_a, size_t number_columns_a,
            size_t number_rows_b, size_t number_columns_b, size_t number_rows_c, size_t number_columns_c, double alpha,
            const double* a, size_t leading_dimension_a, const double* b, size_t leading_dimension_b, double beta, double* c,
            size_t leading_dimension_c) {
         const int m = static_cast<int>(number_rows_c);
         const int n = static_cast<int>(number_columns_c);
         const int k = static_cast<int>(transa == 'N' ? number_columns_a : number_rows_a);
         if (transa == 'N') {
            assert(number_rows_c == number_rows_a);
         }
         else {
            assert(number_rows_c == number_columns_a);
         }
         if (transb == 'N') {
            assert(number_columns_c == number_columns_b);
         }
         else {
            assert(number_columns_c == number_rows_b);
         }
         const int lda = static_cast<int>(leading_dimension_a);
         const int ldb = static_cast<int>(leading_dimension_b);
         const int ldc = static_cast<int>(leading_dimension_c);
         assert(lda >= std::max(1, transa == 'N' ? m : k));
         assert(ldb >= std::max(1, transb == 'N' ? k : n));
         assert(ldc >= std::max(1, m));
         dgemm(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
      }

      // performs back-solve with lower-triangular A:
      // op(A) X = alpha B    or
      // X op(A) = alpha B
      // where
      // op(A) = A   or   op(A) = A^T
      // X is overwritten on B
      inline void triangular_back_solve(char side, char uplo, char transa, char diag, size_t number_rows_b, size_t number_columns_b,
            double alpha, const double* a, size_t leading_dimension_a, double* b, size_t leading_dimension_b) {
         const int m = static_cast<int>(number_rows_b);
         const int n = static_cast<int>(number_columns_b);
         const int lda = static_cast<int>(leading_dimension_a);
         const int ldb = static_cast<int>(leading_dimension_b);
         assert(lda >= std::max(1, side == 'L' ? m : n));
         assert(ldb >= std::max(1, m));
         dtrsm(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
      }

      /*
      // performs one of the matrix-matrix operations
      // B := alpha op(A) B   or   B := alpha B op(A)
      // where
      // alpha is a scalar, B is an m by n matrix, A is a unit, or non-unit, triangular matrix, and
      // op(X) = X   or   op(X) = X^T
      inline void BLAS_triangular_matrix_matrix_product(char side, char uplo, char transa, char diag, size_t number_rows_b,
            size_t number_columns_b, double alpha, double* a, size_t leading_dimension_a, double* b, size_t leading_dimension_b) {
         const int m = static_cast<int>(number_rows_b);
         const int n = static_cast<int>(number_columns_b);
         const int lda = static_cast<int>(leading_dimension_a);
         const int ldb = static_cast<int>(leading_dimension_b);
         assert(lda >= std::max(1, side == 'L' ? m : n));
         assert(ldb >= std::max(1, m));
         dtrmm(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
      }
      */

      // performs symmetric rank-1 update:
      // A := alpha x x^T + A
      inline void symmetric_rank_1_update(char uplo, size_t dimension, double alpha, const double* x, double* a,
            size_t leading_dimension_a) {
         const int n = static_cast<int>(dimension);
         constexpr int increment = 1;
         const int lda = static_cast<int>(leading_dimension_a);
         assert(lda >= std::max(1, n));
         dsyr(&uplo, &n, &alpha, x, &increment, a, &lda);
      }

      // performs symmetric rank k update:
      // C = alpha A A^T + beta C    or
      // C = alpha A^T A + beta C
      inline void symmetric_rank_k_update(char uplo, char trans, size_t dimension_c, size_t number_rows_a,
            size_t number_columns_a, double alpha, const double* a, size_t leading_dimension_a, double beta, double* c,
            size_t leading_dimension_c) {
         const int n = static_cast<int>(dimension_c);
         const int k = static_cast<int>(trans == 'N' ? number_columns_a : number_rows_a);
         const int lda = static_cast<int>(leading_dimension_a);
         assert(lda >= std::max(1, trans == 'N' ? n : k));
         const int ldc = static_cast<int>(leading_dimension_c);
         dsyrk(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
      }
   }
} // namespace

#endif // UNO_BLAS_H