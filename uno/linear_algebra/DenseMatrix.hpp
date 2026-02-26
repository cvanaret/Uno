// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DENSEMATRIX_H
#define UNO_DENSEMATRIX_H

#include <cassert>
#include <ostream>
#include <vector>
#include "linear_algebra/BLAS.hpp"
#include "VectorView.hpp"
#include "symbolic/Inverse.hpp"
#include "symbolic/Multiplication.hpp"
#include "symbolic/Range.hpp"
#include "symbolic/Sum.hpp"
#include "symbolic/Transpose.hpp"

namespace uno {
   // DenseMatrix is an m x n matrix in column-major order where the columns are concatenated in a long vector
   template <typename T>
   class DenseMatrix {
   public:
      using value_type = T;

      DenseMatrix(size_t number_rows, size_t number_columns);
      DenseMatrix& operator=(const DenseMatrix& other);
      ~DenseMatrix() = default;

      // specialized operator= for C := beta * C + A*B^T
      template <typename Matrix1, typename Matrix2, typename Matrix3>
      DenseMatrix& operator=(Sum<ScalarMultiple<Matrix1>, Multiplication<Matrix2, Transpose<Matrix3>>>&& expression);

      // specialized operator= for B := B A⁻ᵀ
      template <typename Matrix>
      DenseMatrix& operator*=(Transpose<Inverse<Matrix>>&& expression);

      [[nodiscard]] T& entry(size_t row_index, size_t column_index);
      [[nodiscard]] const T& entry(size_t row_index, size_t column_index) const;
      // vector view
      [[nodiscard]] MutableVectorView<std::vector<double>> column(size_t column_index);
      [[nodiscard]] VectorView<std::vector<double>> column(size_t column_index) const;

      [[nodiscard]] T* data();
      [[nodiscard]] const T* data() const;
      void fill(T value);
      void print(std::ostream& stream) const;

      const size_t number_rows, number_columns;

   protected:
      std::vector<T> matrix; // column-major ordering
   };

   template <typename T>
   DenseMatrix<T>::DenseMatrix(size_t number_rows, size_t number_columns):
         number_rows(number_rows), number_columns(number_columns),
         matrix(number_rows * number_columns, T(0)) {
   }

   // copy an existing matrix into this object
   template <typename T>
   DenseMatrix<T>& DenseMatrix<T>::operator=(const DenseMatrix& other) {
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
   DenseMatrix<T>& DenseMatrix<T>::operator=(Sum<ScalarMultiple<Matrix1>, Multiplication<Matrix2,
         Transpose<Matrix3>>>&& expression) {
      const double beta = expression.get_left().get_factor();
      const auto& C = expression.get_left().get_expression();
      const auto& A = expression.get_right().get_left();
      const auto& B = expression.get_right().get_right().get_matrix();
      // check that "this" is indeed on the left side of the sum
      if (&C != this) {
         throw std::runtime_error("DenseMatrix::operator= is not defined when C != this");
      }
      if (this->number_rows != A.number_rows || this->number_columns != B.number_rows) {
         throw std::runtime_error("Dimension mismatch in DenseMatrix::operator=");
      }
      const char transa = 'N'; // A
      const char transb = 'T'; // B^T
      const int m = static_cast<int>(A.number_rows); // number of rows of A
      const int n = static_cast<int>(B.number_rows); // number of columns of B^T/number of rows of B
      const int k = static_cast<int>(A.number_columns); // number of columns of A
      constexpr double alpha = 1.;
      const int lda = static_cast<int>(A.number_rows); // leading dimension of A
      const int ldb = static_cast<int>(B.number_rows); // leading dimension of B
      const int ldc = static_cast<int>(this->number_rows); // leading dimension of C/this
      BLAS_matrix_matrix_product(&transa, &transb, &m, &n, &k, &alpha, A.data(), &lda, B.data(), &ldb, &beta,
         this->data(), &ldc);
      return *this;
   }

   // specialized operator= for B *= A⁻ᵀ (solve X Aᵀ := B and overwrite B with X)
   template <typename T>
   template <typename Matrix>
   DenseMatrix<T>& DenseMatrix<T>::operator*=(Transpose<Inverse<Matrix>>&& expression) {
      const auto& A = expression.get_matrix().get_matrix();
      if (this->number_columns != A.number_columns) {
         throw std::runtime_error("Dimension mismatch in DenseMatrix::operator=");
      }
      char transa = 'T'; // op(A) = A^T
      char side = 'R'; // X A^T = alpha B
      char uplo = 'L'; // A is lower triangular
      char diag = 'N';
      int m = static_cast<int>(this->number_rows); // number of rows of B/this
      int n = static_cast<int>(this->number_columns); // number of columns of B/this
      double alpha = 1.;
      int lda = static_cast<int>(A.number_rows); // leading dimension of A
      int ldb = static_cast<int>(this->number_rows); // leading dimension of B/this
      BLAS_triangular_back_solve(&side, &uplo, &transa, &diag, &m, &n, &alpha, A.data(), &lda, this->data(), &ldb);
      return *this;
   }

   template <typename T>
   T& DenseMatrix<T>::entry(size_t row_index, size_t column_index) {
      return this->matrix[column_index * this->number_rows + row_index];
   }

   template <typename T>
   const T& DenseMatrix<T>::entry(size_t row_index, size_t column_index) const {
      return this->matrix[column_index * this->number_rows + row_index];
   }

   template <typename T>
   MutableVectorView<std::vector<double>>  DenseMatrix<T>::column(size_t column_index) {
      return {this->matrix, column_index * this->number_rows, (column_index + 1) * this->number_rows};
   }

   template <typename T>
   VectorView<std::vector<double>>  DenseMatrix<T>::column(size_t column_index) const {
      return {this->matrix, column_index * this->number_rows, (column_index + 1) * this->number_rows};
   }

   template <typename T>
   T* DenseMatrix<T>::data() {
      return this->matrix.data();
   }

   template <typename T>
   const T* DenseMatrix<T>::data() const {
      return this->matrix.data();
   }

   template <typename T>
   void DenseMatrix<T>::fill(T value) {
      for (size_t column_index: Range(this->number_columns)) {
         for (size_t row_index: Range(this->number_rows)) {
            this->entry(row_index, column_index) = value;
         }
      }
   }

   template <typename T>
   void DenseMatrix<T>::print(std::ostream& stream) const {
      stream << "Dense matrix (" << this->number_rows << "x" << this->number_columns << ")\n";
      for (size_t column_index: Range(this->number_columns)) {
         stream << "Column " << column_index << ":";
         for (size_t row_index: Range(this->number_rows)) {
            stream << ' ' << this->entry(row_index, column_index);
         }
         stream << '\n';
      }
   }

   template <typename T>
   std::ostream& operator<<(std::ostream& stream, const DenseMatrix<T>& matrix) {
      matrix.print(stream);
      return stream;
   }
} // namespace

#endif // UNO_DENSEMATRIX_H