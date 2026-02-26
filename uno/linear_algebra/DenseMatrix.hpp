// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DENSEMATRIX_H
#define UNO_DENSEMATRIX_H

#include <ostream>
#include <vector>
#include "linear_algebra/BLASMatrix.hpp"
#include "VectorView.hpp"
#include "symbolic/Range.hpp"

namespace uno {
   // DenseMatrix is an m x n matrix in column-major order where the columns are concatenated in a long vector
   template <typename T>
   class DenseMatrix: public BLASMatrix<T> {
   public:
      using value_type = T;

      // define a Submatrix class: a view of the current matrix with fewer rows and columns (starting from 0)
      // this is for automatizing the calls to BLAS in which dimensions and leading dimension must be provided
      class Submatrix: public BLASMatrix<T> {
      public:
         using value_type = T;

         Submatrix(DenseMatrix<T>& matrix, size_t number_rows, size_t number_columns):
               // the leading dimension is the number of rows of the original matrix
               BLASMatrix<T>(number_rows, number_columns, matrix.number_rows),
               matrix(matrix) {
         }

         // operators
         using BLASMatrix<T>::operator=;
         using BLASMatrix<T>::operator*=;
         // delegate copy operator to that of BLASMatrix
         Submatrix& operator=(const Submatrix& other) {
            BLASMatrix<T>::operator=(other);
            return *this;
         }

         [[nodiscard]] T* data() override {
            return this->matrix.data();
         }

         [[nodiscard]] const T* data() const override {
            return this->matrix.data();
         }

         DenseMatrix<T>& matrix;
      };

      DenseMatrix(size_t number_rows, size_t number_columns);
      ~DenseMatrix() override = default;

      // operators
      using BLASMatrix<T>::operator=;
      using BLASMatrix<T>::operator*=;

      [[nodiscard]] T& entry(size_t row_index, size_t column_index);
      [[nodiscard]] const T& entry(size_t row_index, size_t column_index) const;
      // vector view
      [[nodiscard]] MutableVectorView<std::vector<double>> column(size_t column_index);
      [[nodiscard]] VectorView<std::vector<double>> column(size_t column_index) const;
      [[nodiscard]] Submatrix submatrix(size_t number_rows, size_t number_columns);

      [[nodiscard]] T* data() override;
      [[nodiscard]] const T* data() const override;
      void fill(T value);
      void print(std::ostream& stream) const;

   protected:
      std::vector<T> matrix; // column-major ordering
   };

   template <typename T>
   DenseMatrix<T>::DenseMatrix(size_t number_rows, size_t number_columns):
         BLASMatrix<T>(number_rows, number_columns, number_rows),
         matrix(number_rows * number_columns, T(0)) {
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
   MutableVectorView<std::vector<double>> DenseMatrix<T>::column(size_t column_index) {
      return {this->matrix, column_index * this->number_rows, (column_index + 1) * this->number_rows};
   }

   template <typename T>
   VectorView<std::vector<double>> DenseMatrix<T>::column(size_t column_index) const {
      return {this->matrix, column_index * this->number_rows, (column_index + 1) * this->number_rows};
   }

   template <typename T>
   typename DenseMatrix<T>::Submatrix DenseMatrix<T>::submatrix(size_t number_rows, size_t number_columns) {
      return {*this, number_rows, number_columns};
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