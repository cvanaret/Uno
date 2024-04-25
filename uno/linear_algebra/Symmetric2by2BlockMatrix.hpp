// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRIC2BY2BLOCKMATRIX_H
#define UNO_SYMMETRIC2BY2BLOCKMATRIX_H

#include "VectorView.hpp"
#include "expression/TransposedMatrix.hpp"
#include "linear_algebra/Vector.hpp"

// symmetric block matrix
// (let's start with 2 rows and columns)
// (A    B)
// (B^T  C)
// represented in upper triangular form as the tuple (A, B, C)
template <typename TopLeftBlock, typename TopRightBlock, typename BottomRightBlock>
class Symmetric2by2BlockMatrix {
public:
   Symmetric2by2BlockMatrix(TopLeftBlock&& A, TopRightBlock&& B, BottomRightBlock&& C);
   void for_each(const std::function<void(size_t, size_t, double)>& f) const;
   void product(const std::vector<double>& vector, std::vector<double>& result) const;

protected:
   TopLeftBlock A;
   TopRightBlock B;
   BottomRightBlock C;
};

template <typename TopLeftBlock, typename TopRightBlock, typename BottomRightBlock>
Symmetric2by2BlockMatrix<TopLeftBlock, TopRightBlock, BottomRightBlock>::Symmetric2by2BlockMatrix(TopLeftBlock&& A, TopRightBlock&& B,
      BottomRightBlock&& C): A(std::forward<TopLeftBlock>(A)), B(std::forward<TopRightBlock>(B)), C(std::forward<BottomRightBlock>(C)) {
   // check dimensions
   if (this->A.number_rows() != this->B.number_rows() || this->B.number_columns() != this->C.number_columns()) {
      throw std::runtime_error("The symmetric block matrix was built with blocks of inconsistent sizes.");
   }
}

// iterator
template <typename TopLeftBlock, typename TopRightBlock, typename BottomRightBlock>
void Symmetric2by2BlockMatrix<TopLeftBlock, TopRightBlock, BottomRightBlock>::for_each(const std::function<void(size_t, size_t, double)>& f) const {
   const size_t column_offset = this->A.number_columns();
   const size_t row_offset = this->A.number_rows();

   // top left block
   this->A.for_each(f);

   // top right block
   this->B.for_each([=](size_t row_index, size_t column_index, double entry) {
      f(row_index, column_index + column_offset, entry);
   });

   // bottom left block
   transpose(this->B).for_each([=](size_t row_index, size_t column_index, double entry) {
      f(row_index + row_offset, column_index, entry);
   });

   // bottom right block
   this->C.for_each([=](size_t row_index, size_t column_index, double entry) {
      f(row_index + row_offset, column_index + column_offset, entry);
   });
}

template <typename TopLeftBlock, typename TopRightBlock, typename BottomRightBlock>
void Symmetric2by2BlockMatrix<TopLeftBlock, TopRightBlock, BottomRightBlock>::product(const std::vector<double>& vector, std::vector<double>& result) const {
   initialize_vector(result, 0.);
   // create vector views of the vector and the result
   const VectorView vector_top_part = view(vector, 0, this->A.number_columns());
   const VectorView vector_bottom_part = view(vector, this->A.number_rows(), result.size());
   VectorView result_top_part = view(result, 0, this->A.number_rows());
   VectorView result_bottom_part = view(result, this->A.number_rows(), result.size());

   // carry out the matrix-vector products block-wise
   this->A.product(vector_top_part, result_top_part);
   this->B.product(vector_bottom_part, result_top_part);
   transpose(this->B).product(vector_top_part, result_bottom_part);
   this->C.product(vector_bottom_part, result_bottom_part);
}

#endif // UNO_SYMMETRIC2BY2BLOCKMATRIX_H