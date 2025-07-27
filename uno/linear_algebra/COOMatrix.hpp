// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_COOMATRIX_H
#define UNO_COOMATRIX_H

#include <iostream>
#include <tuple>
#include "symbolic/VectorView.hpp"

namespace uno {
   template <typename IndexType>
   class COOMatrix {
   public:
      IndexType* row_indices;
      IndexType* column_indices;
      double* values;

      COOMatrix(IndexType* row_indices, IndexType* column_indices, double* values):
         row_indices(row_indices), column_indices(column_indices), values(values) { }

      std::tuple<IndexType, IndexType, double> operator[](size_t nonzero_index) const {
         return {this->row_indices[nonzero_index], this->column_indices[nonzero_index], this->values[nonzero_index]};
      }
   };

   template <typename IndexType>
   class COOMatrixView {
   public:
      COOMatrixView(COOMatrix<IndexType>& matrix, size_t start, size_t end):
         matrix(matrix), start(start), end(end) { }

      void print() const {
         std::cout << view(this->matrix.row_indices, this->start, this->end) << '\n';
         std::cout << view(this->matrix.column_indices, this->start, this->end) << '\n';
         std::cout << view(this->matrix.values, this->start, this->end) << '\n';
      }

   protected:
      COOMatrix<IndexType> matrix;
      const size_t start, end;
   };

   template <typename IndexType>
   COOMatrixView<IndexType> view(COOMatrix<IndexType>& matrix, size_t start, size_t end) {
      return {matrix, start, end};
   }
} // namespace

#endif // UNO_COOMATRIX_H