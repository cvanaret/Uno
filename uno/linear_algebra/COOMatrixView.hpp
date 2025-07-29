// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_COOMATRIXBLOCKVIEW_H
#define UNO_COOMATRIXBLOCKVIEW_H

#include <iostream>
#include "symbolic/VectorView.hpp"

namespace uno {
   template <typename IndexType>
   class COOMatrix;

   template <typename IndexType>
   class COOMatrixBlockView {
   public:
      COOMatrixBlockView(COOMatrix<IndexType>& matrix, size_t start, size_t end):
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
   COOMatrixBlockView<IndexType> view(COOMatrix<IndexType>& matrix, size_t start, size_t end) {
      return {matrix, start, end};
   }
} // namespace

#endif // UNO_COOMATRIXBLOCKVIEW_H