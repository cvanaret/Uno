// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DENSEMATRIX_H
#define UNO_DENSEMATRIX_H

#include <iostream>
#include <vector>
#include "symbolic/Range.hpp"

namespace uno {
   // DenseMatrix is an m x n matrix in column-order order where the columns are concatenated in a long vector

   template <typename ElementType>
   class DenseMatrix {
   public:
      DenseMatrix(size_t number_rows, size_t number_columns);
      DenseMatrix(const DenseMatrix& other);
      DenseMatrix& operator=(const DenseMatrix& other);
      ~DenseMatrix() = default;

      ElementType& entry(size_t row_index, size_t column_index);
      const ElementType& entry(size_t row_index, size_t column_index) const;
      ElementType* data();
      void clear();

      //VectorView<const DenseMatrix&> column(size_t column_index) const;

      void print(std::ostream& stream) const;

   protected:
      const size_t number_rows, number_columns;
      std::vector<ElementType> matrix; // column-major order
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
   DenseMatrix<ElementType>& DenseMatrix<ElementType>::operator=(const DenseMatrix& other) {
      for (size_t column_index: Range(this->number_columns)) {
         for (size_t row_index: Range(this->number_rows)) {
            this->entry(row_index, column_index) = other.entry(row_index, column_index);
         }
      }
      return *this;
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