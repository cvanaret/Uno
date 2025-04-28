// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DENSEMATRIX_H
#define UNO_DENSEMATRIX_H

#include <iostream>
#include <vector>
#include "symbolic/Range.hpp"

namespace uno {
   // forward declaration
   template <typename ElementType>
   class DenseMatrix;

   template <typename ElementType>
   class DenseColumn {
   public:
      using value_type = ElementType;

      DenseColumn(const DenseMatrix<ElementType>& matrix, size_t column_index);
      size_t size() const;
      ElementType& operator[](size_t row_index);
      const ElementType& operator[](size_t row_index) const;

   protected:
      const DenseMatrix<ElementType>& matrix;
      const size_t column_index;
   };

   template<typename ElementType>
   DenseColumn<ElementType>::DenseColumn(const DenseMatrix<ElementType> &matrix, size_t column_index):
      matrix(matrix), column_index(column_index) {
   }

   template<typename ElementType>
   size_t DenseColumn<ElementType>::size() const {
      return this->matrix.get_number_rows();
   }

   template<typename ElementType>
   ElementType& DenseColumn<ElementType>::operator[](size_t row_index) {
      return this->matrix.entry(row_index, this->column_index);
   }

   template<typename ElementType>
   const ElementType& DenseColumn<ElementType>::operator[](size_t row_index) const {
      return this->matrix.entry(row_index, this->column_index);
   }

   // DenseMatrix is an m x n matrix in column-order order where the columns are concatenated in a long vector
   template <typename ElementType>
   class DenseMatrix {
   public:
      DenseMatrix(size_t number_rows, size_t number_columns);
      DenseMatrix() = default;
      DenseMatrix(const DenseMatrix& other);
      DenseMatrix& operator=(const DenseMatrix& other) = default;
      DenseMatrix& operator=(DenseMatrix&& other) noexcept = default;
      ~DenseMatrix() = default;

      size_t get_number_rows() const;
      size_t get_number_columns() const;
      ElementType& entry(size_t row_index, size_t column_index);
      const ElementType& entry(size_t row_index, size_t column_index) const;
      // vector view
      DenseColumn<ElementType> column(size_t column_index) const;
      ElementType* data();
      void clear();

      //VectorView<const DenseMatrix&> column(size_t column_index) const;

      void print(std::ostream& stream) const;

   protected:
      size_t number_rows{}, number_columns{};
      std::vector<ElementType> matrix{}; // column-major ordering
   };

   template <typename ElementType>
   DenseMatrix<ElementType>::DenseMatrix(size_t number_rows, size_t number_columns):
         number_rows(number_rows), number_columns(number_columns),
         matrix(number_rows * number_columns, ElementType(0)) {
   }

   template <typename ElementType>
   DenseMatrix<ElementType>::DenseMatrix(const DenseMatrix& other):
         number_rows(other.number_rows), number_columns(other.number_columns),
         matrix(other.number_rows * other.number_columns, ElementType(0)) {
      for (size_t column_index: Range(this->number_columns)) {
         for (size_t row_index: Range(this->number_rows)) {
            this->entry(row_index, column_index) = other.entry(row_index, column_index);
         }
      }
   }

   template <typename ElementType>
   size_t DenseMatrix<ElementType>::get_number_rows() const {
      return this->number_rows;
   }

   template <typename ElementType>
   size_t DenseMatrix<ElementType>::get_number_columns() const {
      return this->number_columns;
   }

   template <typename ElementType>
   ElementType& DenseMatrix<ElementType>::entry(size_t row_index, size_t column_index) {
      return this->matrix[column_index * this->number_rows + row_index];
   }

   template <typename ElementType>
   const ElementType& DenseMatrix<ElementType>::entry(size_t row_index, size_t column_index) const {
      return this->matrix[column_index * this->number_rows + row_index];
   }

   template <typename ElementType>
   DenseColumn<ElementType> DenseMatrix<ElementType>::column(size_t column_index) const {
      return {*this, column_index};
   }

   template <typename ElementType>
   ElementType* DenseMatrix<ElementType>::data() {
      return this->matrix.data();
   }

   template <typename ElementType>
   void DenseMatrix<ElementType>::clear() {
      for (size_t column_index: Range(this->number_columns)) {
         for (size_t row_index: Range(this->number_rows)) {
            this->entry(row_index, column_index) = ElementType(0);
         }
      }
   }

   /*
   template <typename ElementType>
   VectorView<const DenseMatrix&> DenseMatrix<ElementType>::column(size_t column_index) const {
      return {*this, column_index * this->number_rows, column_index * this->number_rows + this->number_rows};
   }
   */

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