// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DENSEMATRIX_H
#define UNO_DENSEMATRIX_H

#include <cstddef>
#include <vector>

namespace uno {
   template <typename ElementType>
   class DenseMatrix {
   public:
      DenseMatrix(size_t number_rows, size_t number_columns);
      DenseMatrix() = default;
      DenseMatrix<ElementType>& operator=(const DenseMatrix<ElementType>& other) = default;
      DenseMatrix<ElementType>& operator=(DenseMatrix<ElementType>&& other) = default;
      ~DenseMatrix() = default;

      ElementType& get(size_t row_index, size_t column_index);
      ElementType* data() const;

   protected:
      const size_t number_rows{}, number_columns{};
      std::vector<ElementType> matrix{}; // column-major ordering
   };

   template <typename ElementType>
   DenseMatrix<ElementType>::DenseMatrix(size_t number_rows, size_t number_columns):
         number_rows(number_rows), number_columns(number_columns),
         matrix(number_rows * number_columns, ElementType(0)) {
   }

   template <typename ElementType>
   ElementType& DenseMatrix<ElementType>::get(size_t row_index, size_t column_index) {
      return this->matrix[column_index * this->number_rows + row_index];
   }

   template <typename ElementType>
   ElementType* DenseMatrix<ElementType>::data() const {
      return this->matrix.data();
   }
} // namespace

#endif // UNO_DENSEMATRIX_H