// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MATRIXVIEW_H
#define UNO_MATRIXVIEW_H

namespace uno {
   // span of an arbitrary container: allocation-free view of a certain length
   class MatrixView {
   public:
      MatrixView(SymmetricMatrix<size_t, double>& matrix, size_t row_start, size_t row_end, size_t column_start, size_t column_end):
            matrix(matrix), row_start(row_start), row_end(row_end), column_start(column_start), column_end(column_end) {
         if (row_end < row_start) {
            throw std::runtime_error("The view ends before its starting point.");
         }
         if (column_end < column_start) {
            throw std::runtime_error("The view ends before its starting point.");
         }
      }

      void insert(double term, size_t row_index, size_t column_index) {
         matrix.insert(term, this->row_start + row_index, this->column_start + column_index);
      }

   protected:
      SymmetricMatrix<size_t, double>& matrix;
      const size_t row_start, row_end, column_start, column_end;
   };

   // free function
   MatrixView view(SymmetricMatrix<size_t, double>& matrix, size_t row_start, size_t row_end, size_t column_start, size_t column_end) {
      return {matrix, row_start, row_end, column_start, column_end};
   }
} // namespace

#endif //UNO_MATRIXVIEW_H
