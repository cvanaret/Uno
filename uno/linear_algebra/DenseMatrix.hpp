// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DENSEMATRIX_H
#define UNO_DENSEMATRIX_H

#include <cassert>
#include <ostream>
#include <vector>
#include "VectorView.hpp"
#include "symbolic/Range.hpp"

namespace uno {
   enum class MatrixShape {
      GENERAL,
      LOWER_TRIANGULAR,
      UPPER_TRIANGULAR
   };

   // DenseMatrix is an m x n matrix in column-major order where the columns are concatenated in a long vector
   template <typename ElementType, MatrixShape Shape = MatrixShape::GENERAL>
   class DenseMatrix {
   public:
      using value_type = ElementType;

      DenseMatrix(size_t number_rows, size_t number_columns);
      DenseMatrix() = default;
      DenseMatrix(const DenseMatrix& other);
      DenseMatrix& operator=(const DenseMatrix& other) = default;
      DenseMatrix& operator=(DenseMatrix&& other) noexcept = default;
      ~DenseMatrix() = default;

      [[nodiscard]] size_t get_number_rows() const;
      [[nodiscard]] size_t get_number_columns() const;
      [[nodiscard]] ElementType& entry(size_t row_index, size_t column_index);
      [[nodiscard]] const ElementType& entry(size_t row_index, size_t column_index) const;
      // vector view
      [[nodiscard]] MutableVectorView<std::vector<double>> column(size_t column_index);
      [[nodiscard]] VectorView<std::vector<double>> column(size_t column_index) const;

      [[nodiscard]] ElementType* data();
      void fill(ElementType value);
      void print(std::ostream& stream) const;

   protected:
      size_t number_rows{}, number_columns{};
      std::vector<ElementType> matrix{}; // column-major ordering
   };

   template <typename ElementType, MatrixShape Shape>
   DenseMatrix<ElementType, Shape>::DenseMatrix(size_t number_rows, size_t number_columns):
         number_rows(number_rows), number_columns(number_columns),
         matrix(number_rows * number_columns, ElementType(0)) {
      if ((Shape == MatrixShape::LOWER_TRIANGULAR || Shape == MatrixShape::UPPER_TRIANGULAR) && number_rows != number_columns) {
         throw std::runtime_error("DenseMatrix: a triangular matrix must be square");
      }
   }

   template <typename ElementType, MatrixShape Shape>
   DenseMatrix<ElementType, Shape>::DenseMatrix(const DenseMatrix& other):
         number_rows(other.number_rows), number_columns(other.number_columns),
         matrix(other.number_rows * other.number_columns, ElementType(0)) {
      for (size_t column_index: Range(this->number_columns)) {
         for (size_t row_index: Range(this->number_rows)) {
            this->entry(row_index, column_index) = other.entry(row_index, column_index);
         }
      }
   }

   template <typename ElementType, MatrixShape Shape>
   size_t DenseMatrix<ElementType, Shape>::get_number_rows() const {
      return this->number_rows;
   }

   template <typename ElementType, MatrixShape Shape>
   size_t DenseMatrix<ElementType, Shape>::get_number_columns() const {
      return this->number_columns;
   }

   template <typename ElementType, MatrixShape Shape>
   ElementType& DenseMatrix<ElementType, Shape>::entry(size_t row_index, size_t column_index) {
      assert(row_index < this->get_number_rows() && column_index < this->get_number_columns() &&
         "DenseMatrix::entry: indices out of bounds");
      return this->matrix[column_index * this->number_rows + row_index];
   }

   template <typename ElementType, MatrixShape Shape>
   const ElementType& DenseMatrix<ElementType, Shape>::entry(size_t row_index, size_t column_index) const {
      return this->matrix[column_index * this->number_rows + row_index];
   }

   template <typename ElementType, MatrixShape Shape>
   MutableVectorView<std::vector<double>>  DenseMatrix<ElementType, Shape>::column(size_t column_index) {
      return {this->matrix, column_index * this->number_rows, (column_index + 1) * this->number_rows};
   }

   template <typename ElementType, MatrixShape Shape>
   VectorView<std::vector<double>>  DenseMatrix<ElementType, Shape>::column(size_t column_index) const {
      return {this->matrix, column_index * this->number_rows, (column_index + 1) * this->number_rows};
   }

   template <typename ElementType, MatrixShape Shape>
   ElementType* DenseMatrix<ElementType, Shape>::data() {
      return this->matrix.data();
   }

   template <typename ElementType, MatrixShape Shape>
   void DenseMatrix<ElementType, Shape>::fill(ElementType value) {
      for (size_t column_index: Range(this->number_columns)) {
         for (size_t row_index: Range(this->number_rows)) {
            this->entry(row_index, column_index) = value;
         }
      }
   }

   template <typename ElementType, MatrixShape Shape>
   void DenseMatrix<ElementType, Shape>::print(std::ostream& stream) const {
      stream << "Dense matrix (" << this->number_rows << "x" << this->number_columns << ")\n";
      for (size_t column_index: Range(this->number_columns)) {
         stream << "Column " << column_index << ":";
         for (size_t row_index: Range(this->number_rows)) {
            stream << ' ' << this->entry(row_index, column_index);
         }
         stream << '\n';
      }
   }

   template <typename ElementType, MatrixShape Shape>
   std::ostream& operator<<(std::ostream& stream, const DenseMatrix<ElementType, Shape>& matrix) {
      matrix.print(stream);
      return stream;
   }
} // namespace

#endif // UNO_DENSEMATRIX_H