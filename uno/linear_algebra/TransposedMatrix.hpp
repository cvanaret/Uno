// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_TRANSPOSEDMATRIX_H
#define UNO_TRANSPOSEDMATRIX_H

#include "Matrix.hpp"

template <typename MatrixType>
class TransposedMatrix: Matrix<typename std::remove_reference<MatrixType>::type::value_type> {
public:
   explicit TransposedMatrix(MatrixType&& matrix);
   void for_each(const std::function<void(size_t, size_t, typename TransposedMatrix::value_type)>& f) const;
   template <typename VectorType, typename ResultType>
   void product(const VectorType& vector, ResultType& result) const;

protected:
   MatrixType matrix;
};

template <typename MatrixType>
TransposedMatrix<MatrixType>::TransposedMatrix(MatrixType&& matrix): matrix(std::forward<MatrixType>(matrix)) {
}

template <typename MatrixType>
void TransposedMatrix<MatrixType>::for_each(const std::function<void(size_t, size_t, typename TransposedMatrix::value_type)>& f) const {
   this->matrix.for_each([&](size_t row_index, size_t column_index, typename TransposedMatrix::value_type entry) {
      // switch row and column indices
      f(column_index, row_index, entry);
   });
}

template <typename MatrixType>
template <typename VectorType, typename ResultType>
void TransposedMatrix<MatrixType>::product(const VectorType& vector, ResultType& result) const {
   this->for_each([&](size_t row_index, size_t column_index, double entry) {
      result[row_index] += entry*vector[column_index];
      // off-diagonal term on the other side of the diagonal
      if (row_index != column_index) {
         result[column_index] += entry*vector[row_index];
      }
   });
}

// free function
template <typename MatrixType>
TransposedMatrix<MatrixType> transpose(MatrixType&& matrix) {
   return TransposedMatrix(std::forward<MatrixType>(matrix));
}

#endif // UNO_TRANSPOSEDMATRIX_H