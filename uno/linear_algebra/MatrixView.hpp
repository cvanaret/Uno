// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MATRIXVIEW_H
#define UNO_MATRIXVIEW_H

#include <cassert>
#include "Matrix.hpp"

namespace uno {
   // abstract class
   template <typename IndexType, typename ElementType>
   class MatrixView: public Matrix<IndexType, ElementType> {
   public:
      MatrixView(Matrix<IndexType, ElementType>& source_matrix, IndexType start_row, IndexType start_column, IndexType end_row, IndexType end_column,
         size_t number_nonzeros);
      virtual ~MatrixView() = default;

      void insert(ElementType term, IndexType row_index, IndexType column_index) override;
      void finalize_column(IndexType column_index) override;

   protected:
      Matrix<IndexType, ElementType>& source_matrix;
      const IndexType start_row, start_column, end_row, end_column;
      const size_t number_nonzeros;
      size_t current_number_nonzeros{0};
   };

   template <typename IndexType, typename ElementType>
   MatrixView<IndexType, ElementType>::MatrixView(Matrix<IndexType, ElementType> &source_matrix, IndexType start_row, IndexType start_column,
         IndexType end_row, IndexType end_column, size_t number_nonzeros):
         Matrix<IndexType, ElementType>(end_row - start_row, end_column - start_column, number_nonzeros),
         source_matrix(source_matrix), start_row(start_row), start_column(start_column), end_row(end_row),
         end_column(end_column), number_nonzeros(number_nonzeros) {
      assert(this->number_nonzeros <= this->source_matrix.get_number_nonzeros() && "The capacity of the matrix view exceeds that of the source matrix");
   }

   template <typename IndexType, typename ElementType>
   void MatrixView<IndexType, ElementType>::insert(ElementType term, IndexType row_index, IndexType column_index) {
      assert(this->current_number_nonzeros < this->number_nonzeros && "The capacity of the matrix view has been reached");
      assert(row_index < this->end_row && "The row index exceeds the size of the matrix view");
      assert(column_index < this->end_column && "The column index exceeds the size of the matrix view");

      this->source_matrix.insert(term, this->start_row + row_index, this->start_column + column_index);
      this->current_number_nonzeros++;
   }

   template <typename IndexType, typename ElementType>
   void MatrixView<IndexType, ElementType>::finalize_column(IndexType column_index) {
      this->source_matrix.finalize_column(this->start_column + column_index);
   }
} // namespace

#endif // UNO_MATRIXVIEW_H