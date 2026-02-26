// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BLASMATRIX_H
#define UNO_BLASMATRIX_H

#include <cassert>
#include "linear_algebra/BLAS.hpp"
#include "VectorView.hpp"
#include "symbolic/Inverse.hpp"
#include "symbolic/Multiplication.hpp"
#include "symbolic/Sum.hpp"
#include "symbolic/Transpose.hpp"

namespace uno {
   template <typename T>
   class BLASMatrix {
   public:
      using value_type = T;

      BLASMatrix(size_t number_rows, size_t number_columns, size_t leading_dimension);
      virtual ~BLASMatrix() = default;

      // copy an existing matrix into this object
      BLASMatrix<T>& operator=(const BLASMatrix& other);

      // specialized operator= for C := beta * C + A*B^T
      template <typename Matrix1, typename Matrix2, typename Matrix3>
      BLASMatrix& operator=(Sum<ScalarMultiple<Matrix1>, Multiplication<Matrix2, Transpose<Matrix3>>>&& expression);

      // specialized operator= for B := B A⁻ᵀ
      template <typename Matrix>
      BLASMatrix& operator*=(Transpose<Inverse<Matrix>>&& expression);

      [[nodiscard]] virtual T* data() = 0;
      [[nodiscard]] virtual const T* data() const = 0;

      const size_t number_rows, number_columns, leading_dimension;
   };

   template <typename T>
   BLASMatrix<T>::BLASMatrix(size_t number_rows, size_t number_columns, size_t leading_dimension):
         number_rows(number_rows), number_columns(number_columns), leading_dimension(leading_dimension) {
   }

   // copy an existing matrix into this object
   template <typename T>
   BLASMatrix<T>& BLASMatrix<T>::operator=(const BLASMatrix& other) {
      assert(other.number_rows == this->number_rows && other.number_columns == this->number_columns);
      const int size = static_cast<int>(this->number_rows * this->number_columns);
      constexpr int increment = 1;
      BLAS_copy_vector(&size, other.data(), &increment, this->data(), &increment);
      return *this;
   }

   // specialized operator= for C := beta * C + A*B^T
   // use different matrix types in case one of them has a different type (e.g., submatrix)
   template <typename T>
   template <typename Matrix1, typename Matrix2, typename Matrix3>
   BLASMatrix<T>& BLASMatrix<T>::operator=(Sum<ScalarMultiple<Matrix1>, Multiplication<Matrix2,
         Transpose<Matrix3>>>&& expression) {
      const double beta = expression.get_left().get_factor();
      const auto& C = expression.get_left().get_expression();
      const auto& A = expression.get_right().get_left();
      const auto& B = expression.get_right().get_right().get_matrix();
      // check that "this" is indeed on the left side of the sum
      if (&C != this) {
         throw std::runtime_error("BLASMatrix::operator= is not defined when C != this");
      }
      if (this->number_rows != A.number_rows || this->number_columns != B.number_rows) {
         throw std::runtime_error("Dimension mismatch in BLASMatrix::operator=");
      }
      const char transa = 'N'; // A
      const char transb = 'T'; // B^T
      const int m = static_cast<int>(A.number_rows); // number of rows of A
      const int n = static_cast<int>(B.number_rows); // number of columns of B^T/number of rows of B
      const int k = static_cast<int>(A.number_columns); // number of columns of A
      constexpr double alpha = 1.;
      const int lda = static_cast<int>(A.leading_dimension); // leading dimension of A
      const int ldb = static_cast<int>(B.leading_dimension); // leading dimension of B
      const int ldc = static_cast<int>(this->number_rows); // leading dimension of C/this
      BLAS_matrix_matrix_product(&transa, &transb, &m, &n, &k, &alpha, A.data(), &lda, B.data(), &ldb, &beta,
         this->data(), &ldc);
      return *this;
   }

   // specialized operator= for B *= A⁻ᵀ (solve X Aᵀ := B and overwrite B with X)
   template <typename T>
   template <typename Matrix>
   BLASMatrix<T>& BLASMatrix<T>::operator*=(Transpose<Inverse<Matrix>>&& expression) {
      const auto& A = expression.get_matrix().get_matrix();
      if (this->number_columns != A.number_columns) {
         throw std::runtime_error("Dimension mismatch in BLASMatrix::operator*=");
      }
      char transa = 'T'; // op(A) = A^T
      char side = 'R'; // X A^T = alpha B
      char uplo = 'L'; // A is lower triangular
      char diag = 'N';
      int m = static_cast<int>(this->number_rows); // number of rows of B/this
      int n = static_cast<int>(this->number_columns); // number of columns of B/this
      double alpha = 1.;
      int lda = static_cast<int>(A.leading_dimension); // leading dimension of A
      int ldb = static_cast<int>(this->leading_dimension); // leading dimension of B/this
      BLAS_triangular_back_solve(&side, &uplo, &transa, &diag, &m, &n, &alpha, A.data(), &lda, this->data(), &ldb);
      return *this;
   }
} // namespace

#endif // UNO_BLASMATRIX_H