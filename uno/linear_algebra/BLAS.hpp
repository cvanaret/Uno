// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BLAS_H
#define UNO_BLAS_H

#include "fortran_interface.h"
#define dcopy FC_GLOBAL_(dcopy, DCOPY)
#define dscal FC_GLOBAL_(dscal, DSCAL)
#define daxpy FC_GLOBAL_(daxpy, DAXPY)
#define ddot FC_GLOBAL_(ddot, DDOT)
#define dgemv FC_GLOBAL_(dgemv, DGEMV)
#define BLAS_triangular_back_solve FC_GLOBAL_(dtrsm, DTRSM)
#define BLAS_matrix_matrix_product FC_GLOBAL_(dgemm, DGEMM)
#define BLAS_triangular_matrix_matrix_product FC_GLOBAL_(dtrmm, DTRMM)

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
   void BLAS_triangular_back_solve(const char* side, const char* uplo, const char* transa, const char* diag, const int* m,
      const int* n, const double* alpha, const double* a, const int* lda, double* b, const int* ldb);

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

namespace uno {
   // y := x
   inline void BLAS_copy_vector(size_t size, const double* x, double* y) {
      const int n = static_cast<int>(size);
      constexpr int increment = 1;
      dcopy(&n, x, &increment, y, &increment);
   }

   // x := alpha x
   inline void BLAS_scale_vector(size_t size, double alpha, double* x) {
      const int n = static_cast<int>(size);
      constexpr int increment = 1;
      dscal(&n, &alpha, x, &increment);
   }

   // y := alpha x + y
   inline void BLAS_add_vectors(size_t size, double alpha, const double* x, double* y) {
      const int n = static_cast<int>(size);
      constexpr int increment = 1;
      daxpy(&n, &alpha, x, &increment, y, &increment);
   }

   // x^T y
   inline double BLAS_dot_product(size_t size, const double* x, const double* y) {
      const int n = static_cast<int>(size);
      constexpr int increment = 1;
      return ddot(&n, x, &increment, y, &increment);
   }

   // performs the matrix-vector operations
   // y := alpha A x + beta y,   or   y := alpha A^T x + beta y,
   inline void BLAS_matrix_vector_product(char trans, size_t number_rows, size_t number_colums, double alpha, const double* a,
         size_t leading_dimension, const double* x, double beta, double* y) {
      const int m = static_cast<int>(number_rows);
      const int n = static_cast<int>(number_colums);
      const int lda = static_cast<int>(leading_dimension);
      constexpr int increment = 1;
      dgemv(&trans, &m, &n, &alpha, a, &lda, x, &increment, &beta, y, &increment);
   }
}

#endif // UNO_BLAS_H