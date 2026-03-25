// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BLASMATRIX_H
#define UNO_BLASMATRIX_H

#include <cassert>
#include "Vector.hpp"
#include "linear_algebra/BLAS.hpp"
#include "linear_algebra/LAPACK.hpp"
#include "VectorView.hpp"
#include "symbolic/Inverse.hpp"
#include "symbolic/Multiplication.hpp"
#include "symbolic/Sum.hpp"
#include "symbolic/Transpose.hpp"
#include "tools/Logger.hpp"

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

      // specialized operator= for C += A^T*B
      template <typename Matrix1, typename Matrix2>
      BLASMatrix<T>& operator+=(Multiplication<Transpose<Matrix1>, Matrix2>&& expression);

      // specialized operator+= for low-rank update C := A A^T
      template <typename Matrix>
      BLASMatrix& operator=(Multiplication<Matrix, Transpose<Matrix>>&& expression);

      // specialized operator+= for low-rank update C += alpha A^T A
      template <typename Matrix>
      BLASMatrix& operator+=(ScalarMultiple<Multiplication<Transpose<Matrix>, Matrix>>&& expression);

      // specialized operator= for B := B A⁻ᵀ
      template <typename Matrix>
      BLASMatrix& operator*=(Transpose<Inverse<Matrix>>&& expression);

      [[nodiscard]] bool compute_cholesky_factorization();
      [[nodiscard]] std::pair<bool, std::vector<int>> compute_bunch_kaufman_factorization();

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
      blas1::copy(this->number_rows * this->number_columns, other.data(), this->data());
      return *this;
   }

   // specialized operator= for C := beta * C + A*B^T
   // use different matrix types in case one of them has a different type (e.g., submatrix)
   template <typename T>
   template <typename Matrix1, typename Matrix2, typename Matrix3>
   BLASMatrix<T>& BLASMatrix<T>::operator=(Sum<ScalarMultiple<Matrix1>, Multiplication<Matrix2,
         Transpose<Matrix3>>>&& expression) {
      const T beta = expression.get_left().get_factor();
      const auto& C = expression.get_left().get_expression();
      const auto& A = expression.get_right().get_left();
      const auto& B = expression.get_right().get_right().get_matrix();
      // check that "this" is indeed on the left side of the sum
      if (&C != this) {
         throw std::runtime_error("BLASMatrix::operator= is not defined when C != this");
      }
      if (this->number_rows != A.number_rows || this->number_columns != B.number_rows) {
         throw std::runtime_error("Dimension mismatch");
      }
      constexpr char transa = 'N'; // A
      constexpr char transb = 'T'; // B^T
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

   // specialized operator= for C += A^T*B
   // use different matrix types in case one of them has a different type (e.g., submatrix)
   template <typename T>
   template <typename Matrix1, typename Matrix2>
   BLASMatrix<T>& BLASMatrix<T>::operator+=(Multiplication<Transpose<Matrix1>, Matrix2>&& expression) {
      const auto& A = expression.get_left().get_matrix();
      const auto& B = expression.get_right();
      if (A.number_rows != B.number_rows) {
         throw std::runtime_error("Dimension mismatch");
      }
      constexpr char transa = 'T'; // A^T
      constexpr char transb = 'N'; // B
      const int m = static_cast<int>(A.number_columns); // number of rows of A^T/number of columns of A
      const int n = static_cast<int>(B.number_columns); // number of columns of B
      const int k = static_cast<int>(A.number_rows); // number of columns of A^T/number of rows of A
      constexpr double alpha = 1.;
      const int lda = static_cast<int>(A.leading_dimension); // leading dimension of A
      const int ldb = static_cast<int>(B.leading_dimension); // leading dimension of B
      constexpr double beta = 1.;
      const int ldc = static_cast<int>(this->number_rows); // leading dimension of C/this
      BLAS_matrix_matrix_product(&transa, &transb, &m, &n, &k, &alpha, A.data(), &lda, B.data(), &ldb, &beta,
         this->data(), &ldc);
      return *this;
   }

   // specialized operator+= for low-rank update C := A A^T
   template <typename T>
   template <typename Matrix>
   BLASMatrix<T>& BLASMatrix<T>::operator=(Multiplication<Matrix, Transpose<Matrix>>&& expression) {
      assert(this->number_rows == this->number_columns);
      // decode expression as alpha A B^T
      const auto& A = expression.get_left();
      const auto& B = expression.get_right().get_matrix();
      const size_t correction_rank = A.number_columns;
      DEBUG << "Performing rank " << correction_rank << " update\n";
      // check that A and B are the same object
      if (&A != &B) {
         throw std::runtime_error("BLASMatrix::operator+=: low-rank update called on two different correction matrices");
      }
      blas3::symmetric_high_rank_update('L', 'N', this->number_rows, A.number_rows, A.number_columns, 1., A.data(),
         A.leading_dimension, 0., this->data(), this->leading_dimension);
      return *this;
   }

   // specialized operator+= for low-rank update C += alpha A^T A
   template <typename T>
   template <typename Matrix>
   BLASMatrix<T>& BLASMatrix<T>::operator+=(ScalarMultiple<Multiplication<Transpose<Matrix>, Matrix>>&& expression) {
      assert(this->number_rows == this->number_columns);
      // decode expression as alpha A^T B
      const T alpha = expression.get_factor();
      const auto& A = expression.get_expression().get_left().get_matrix();
      const auto& B = expression.get_expression().get_right();
      const size_t correction_rank = A.number_columns;
      DEBUG << "Performing rank " << correction_rank << " update\n";
      // check that A and B are the same object
      if (&A != &B) {
         throw std::runtime_error("BLASMatrix::operator+=: low-rank update called on two different correction matrices");
      }
      blas3::symmetric_high_rank_update('L', 'T', this->number_rows, A.number_rows, A.number_columns, alpha, A.data(),
         A.leading_dimension, 1., this->data(), this->leading_dimension);
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
      constexpr char transa = 'T'; // op(A) = A^T
      constexpr char side = 'R'; // X A^T = alpha B
      constexpr char uplo = 'L'; // A is lower triangular
      constexpr char diag = 'N';
      const int m = static_cast<int>(this->number_rows); // number of rows of B/this
      const int n = static_cast<int>(this->number_columns); // number of columns of B/this
      constexpr double alpha = 1.;
      const int lda = static_cast<int>(A.leading_dimension); // leading dimension of A
      const int ldb = static_cast<int>(this->leading_dimension); // leading dimension of B/this
      BLAS_triangular_back_solve(&side, &uplo, &transa, &diag, &m, &n, &alpha, A.data(), &lda, this->data(), &ldb);
      return *this;
   }

   template <typename T>
   bool BLASMatrix<T>::compute_cholesky_factorization() {
      return lapack::cholesky_factorization('L', this->number_rows, this->data(), this->leading_dimension);
   }

   template <typename T>
   std::pair<bool, std::vector<int>> BLASMatrix<T>::compute_bunch_kaufman_factorization() {
      return lapack::bunch_kaufman_factorization('L', this->number_rows, this->data(), this->leading_dimension);
   }

   template <typename T>
   bool solve_bunch_kaufman(const BLASMatrix<T>& matrix, const Vector<double>& rhs, Vector<double>& result,
         const std::vector<int>& ipiv) {
      // copy the RHS into the result
      result = rhs;
      // solve A X = B
      return lapack::bunch_kaufman_solve('L', matrix.number_rows, matrix.data(), matrix.leading_dimension, ipiv.data(),
         result.data());
   }
} // namespace

#endif // UNO_BLASMATRIX_H