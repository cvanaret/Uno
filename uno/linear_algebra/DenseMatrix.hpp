// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DENSEMATRIX_H
#define UNO_DENSEMATRIX_H

#include <cassert>
#include <ostream>
#include <vector>
#include <iostream>
#include "linear_algebra/BLAS.hpp"
#include "VectorView.hpp"
#include "symbolic/Multiplication.hpp"
#include "symbolic/Range.hpp"
#include "symbolic/Sum.hpp"
#include "symbolic/Transpose.hpp"

namespace uno {
   // DenseMatrix is an m x n matrix in column-major order where the columns are concatenated in a long vector
   template <typename ElementType>
   class DenseMatrix {
   public:
      using value_type = ElementType;

      DenseMatrix(size_t number_rows, size_t number_columns);

      // copy an existing matrix into this object
      DenseMatrix& operator=(DenseMatrix& other) {
         assert(other.number_rows == this->number_rows && other.number_columns == this->number_columns);
         const int size = static_cast<int>(this->number_rows * this->number_columns);
         constexpr int increment = 1;
         BLAS_copy_vector(&size, other.data(), &increment, this->data(), &increment);
         return *this;
      }

      ~DenseMatrix() = default;

      // specialized operator= for C := beta * C + A*B^T
      // use different matrix types in case one of them has a different type (e.g., submatrix)
      template <typename Matrix1, typename Matrix2, typename Matrix3>
      DenseMatrix& operator=(Sum<ScalarMultiple<Matrix1>, Multiplication<Matrix2, Transpose<Matrix3>>>&& expression) {
         const double beta = expression.get_left().get_factor();
         auto& C = expression.get_left().get_expression();
         auto& A = expression.get_right().get_left();
         auto& B = expression.get_right().get_right().get_matrix();
         // check that "this" is indeed on the left side of the sum
         if (&C != this) {
            throw std::runtime_error("DenseMatrix::operator= is not defined when C != this");
         }
         if (this->number_rows != A.number_rows || this->number_columns != B.number_rows) {
            throw std::runtime_error("Dimension mismatch in DenseMatrix::operator=");
         }
         /*
         std::cout << "C/this has size (" << this->number_rows << ", " << this->number_columns << ")\n";
         std::cout << "A has size (" << A.number_rows << ", " << A.number_columns << ")\n";
         std::cout << "B has size (" << B.number_rows << ", " << B.number_columns << ")\n";
         */
         const char transa = 'N'; // A
         const char transb = 'T'; // B^T
         const int m = static_cast<int>(A.number_rows); // number of rows of A
         const int n = static_cast<int>(B.number_rows); // number of columns of B^T
         const int k = static_cast<int>(A.number_columns); // number of columns of A
         constexpr double alpha = 1.;
         const int lda = static_cast<int>(A.number_rows); // leading dimension of A
         const int ldb = static_cast<int>(B.number_rows); // leading dimension of B
         const int ldc = static_cast<int>(this->number_rows); // leading dimension of C/this
         BLAS_matrix_matrix_product(&transa, &transb, &m, &n, &k, &alpha, A.data(), &lda, B.data(), &ldb,
            &beta, this->data(), &ldc);
         return *this;
      }

      [[nodiscard]] ElementType& entry(size_t row_index, size_t column_index);
      [[nodiscard]] const ElementType& entry(size_t row_index, size_t column_index) const;
      // vector view
      [[nodiscard]] MutableVectorView<std::vector<double>> column(size_t column_index);
      [[nodiscard]] VectorView<std::vector<double>> column(size_t column_index) const;

      [[nodiscard]] ElementType* data();
      void fill(ElementType value);
      void print(std::ostream& stream) const;

      const size_t number_rows, number_columns;

   protected:
      std::vector<ElementType> matrix; // column-major ordering
   };

   template <typename ElementType>
   DenseMatrix<ElementType>::DenseMatrix(size_t number_rows, size_t number_columns):
         number_rows(number_rows), number_columns(number_columns),
         matrix(number_rows * number_columns, ElementType(0)) {
   }

   template <typename ElementType>
   ElementType& DenseMatrix<ElementType>::entry(size_t row_index, size_t column_index) {
      assert(row_index < this->number_rows && column_index < this->number_columns &&
         "DenseMatrix::entry: indices out of bounds");
      return this->matrix[column_index * this->number_rows + row_index];
   }

   template <typename ElementType>
   const ElementType& DenseMatrix<ElementType>::entry(size_t row_index, size_t column_index) const {
      return this->matrix[column_index * this->number_rows + row_index];
   }

   template <typename ElementType>
   MutableVectorView<std::vector<double>>  DenseMatrix<ElementType>::column(size_t column_index) {
      return {this->matrix, column_index * this->number_rows, (column_index + 1) * this->number_rows};
   }

   template <typename ElementType>
   VectorView<std::vector<double>>  DenseMatrix<ElementType>::column(size_t column_index) const {
      return {this->matrix, column_index * this->number_rows, (column_index + 1) * this->number_rows};
   }

   template <typename ElementType>
   ElementType* DenseMatrix<ElementType>::data() {
      return this->matrix.data();
   }

   template <typename ElementType>
   void DenseMatrix<ElementType>::fill(ElementType value) {
      for (size_t column_index: Range(this->number_columns)) {
         for (size_t row_index: Range(this->number_rows)) {
            this->entry(row_index, column_index) = value;
         }
      }
   }

   template <typename ElementType>
   void DenseMatrix<ElementType>::print(std::ostream& stream) const {
      stream << "Dense matrix (" << this->number_rows << "x" << this->number_columns << ")\n";
      for (size_t column_index: Range(this->number_columns)) {
         stream << "Column " << column_index << ":";
         for (size_t row_index: Range(this->number_rows)) {
            stream << ' ' << this->entry(row_index, column_index);
         }
         stream << '\n';
      }
   }

   template <typename ElementType>
   std::ostream& operator<<(std::ostream& stream, const DenseMatrix<ElementType>& matrix) {
      matrix.print(stream);
      return stream;
   }
} // namespace

#endif // UNO_DENSEMATRIX_H