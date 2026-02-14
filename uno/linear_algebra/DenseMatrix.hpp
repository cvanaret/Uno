// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DENSEMATRIX_H
#define UNO_DENSEMATRIX_H

#include <cassert>
#include <ostream>
#include <vector>
#include "symbolic/Range.hpp"

namespace uno {
   // MatrixType can be const DenseMatrix (read-only) or DenseMatrix (writable)
   template <typename MatrixType>
   class DenseColumn {
   public:
      using value_type = typename std::remove_reference_t<MatrixType>::value_type;

      DenseColumn(MatrixType& matrix, size_t column_index): matrix(matrix), column_index(column_index) {
      }

      template <typename Expression>
      DenseColumn& operator=(Expression&& expression) {
         static_assert(std::is_same_v<typename std::remove_reference_t<Expression>::value_type, value_type>);
         for (size_t index: Range(this->size())) {
            this->operator[](index) = expression[index];
         }
         return *this;
      }

      template <typename Expression>
      DenseColumn& operator+=(Expression&& expression) {
         static_assert(std::is_same_v<typename std::remove_reference_t<Expression>::value_type, value_type>);
         for (size_t index: Range(this->size())) {
            this->operator[](index) += expression[index];
         }
         return *this;
      }

      [[nodiscard]] size_t size() const {
         return this->matrix.get_number_rows();
      }

      [[nodiscard]] value_type& operator[](size_t row_index) {
         return this->matrix.entry(row_index, this->column_index);
      }

      [[nodiscard]] const value_type& operator[](size_t row_index) const {
         return this->matrix.entry(row_index, this->column_index);
      }

      value_type* data() {
         return this->matrix.data() + this->column_index * this->matrix.get_number_rows();
      }

      const value_type* data() const {
         return this->matrix.data() + this->column_index * this->matrix.get_number_rows();
      }

   protected:
      MatrixType& matrix;
      const size_t column_index;
   };

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
      [[nodiscard]] DenseColumn<DenseMatrix<ElementType, Shape>> column(size_t column_index);
      [[nodiscard]] DenseColumn<const DenseMatrix<ElementType, Shape>> column(size_t column_index) const;

      [[nodiscard]] ElementType* data();
      void clear();
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
   DenseColumn<DenseMatrix<ElementType, Shape>> DenseMatrix<ElementType, Shape>::column(size_t column_index) {
      return {*this, column_index};
   }

   template <typename ElementType, MatrixShape Shape>
   DenseColumn<const DenseMatrix<ElementType, Shape>> DenseMatrix<ElementType, Shape>::column(size_t column_index) const {
      return {*this, column_index};
   }

   template <typename ElementType, MatrixShape Shape>
   ElementType* DenseMatrix<ElementType, Shape>::data() {
      return this->matrix.data();
   }

   template <typename ElementType, MatrixShape Shape>
   void DenseMatrix<ElementType, Shape>::clear() {
      for (size_t column_index: Range(this->number_columns)) {
         for (size_t row_index: Range(this->number_rows)) {
            this->entry(row_index, column_index) = ElementType(0);
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