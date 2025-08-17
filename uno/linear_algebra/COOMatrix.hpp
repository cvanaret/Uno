// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_COOMATRIX_H
#define UNO_COOMATRIX_H

#include <tuple>
#include "Matrix.hpp"

namespace uno {
   template <typename IndexType>
   class COOMatrix: public Matrix<IndexType> {
   public:
      IndexType* row_indices;
      IndexType* column_indices;
      double* values;

      COOMatrix(IndexType* row_indices, IndexType* column_indices, double* values):
         row_indices(row_indices), column_indices(column_indices), values(values) { }

      std::tuple<IndexType, IndexType, double> operator[](size_t nonzero_index) const override {
         return {this->row_indices[nonzero_index], this->column_indices[nonzero_index], this->values[nonzero_index]};
      }
   };
} // namespace

#endif // UNO_COOMATRIX_H