// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRIC2BY2BLOCKMATRIX_H
#define UNO_SYMMETRIC2BY2BLOCKMATRIX_H

#include "VectorView.hpp"
#include "TransposedMatrix.hpp"

// symmetric block matrix
// (let's start with 2 rows and columns)
// (A    B)
// (B^T  C)
// represented in upper triangular form as the tuple (A, B, C)
template <typename TopLeftBlock, typename TopRightBlock, typename BottomRightBlock>
class Symmetric2by2BlockMatrix {
public:
   Symmetric2by2BlockMatrix(const TopLeftBlock& top_left_block, const TopRightBlock& top_right_block, const BottomRightBlock& bottom_right_block);
   void for_each(const std::function<void(size_t, size_t, double)>& f) const;
   void product(const std::vector<double>& vector, std::vector<double>& result) const;

protected:
   const TopLeftBlock& top_left_block;
   const TopRightBlock& top_right_block;
   const BottomRightBlock& bottom_right_block;
};

template <typename TopLeftBlock, typename TopRightBlock, typename BottomRightBlock>
Symmetric2by2BlockMatrix<TopLeftBlock, TopRightBlock, BottomRightBlock>::Symmetric2by2BlockMatrix(const TopLeftBlock& top_left_block,
      const TopRightBlock& top_right_block, const BottomRightBlock& bottom_right_block):
      top_left_block(top_left_block), top_right_block(top_right_block), bottom_right_block(bottom_right_block) {
   // check dimensions
   if (this->top_left_block.number_rows() != this->top_right_block.number_rows() ||
      this->top_right_block.number_columns() != this->bottom_right_block.number_columns()) {
      throw std::runtime_error("The symmetric block matrix was built with blocks of inconsistent sizes.");
   }
}

// iterator: goes through upper triangular terms
template <typename TopLeftBlock, typename TopRightBlock, typename BottomRightBlock>
void Symmetric2by2BlockMatrix<TopLeftBlock, TopRightBlock, BottomRightBlock>::for_each(const std::function<void(size_t, size_t, double)>& f) const {
   // top left block
   this->top_left_block.for_each(f);
   // top right block: offset on column indices
   const size_t column_offset = this->top_left_block.number_columns();
   this->top_right_block.for_each([=](size_t row_index, size_t column_index, double entry) {
      f(row_index, column_index + column_offset, entry);
   });
   // bottom right block: offset on row and column indices
   const size_t row_offset = this->top_left_block.number_rows();
   this->bottom_right_block.for_each([=](size_t row_index, size_t column_index, double entry) {
      f(row_index + row_offset, column_index + column_offset, entry);
   });
}

template <typename TopLeftBlock, typename TopRightBlock, typename BottomRightBlock>
void Symmetric2by2BlockMatrix<TopLeftBlock, TopRightBlock, BottomRightBlock>::product(const std::vector<double>& vector, std::vector<double>& result) const {
   // create vector views (= lazy subvectors) of the vector and the result
   const VectorView vector_top_part = view(vector, 0, this->top_left_block.number_columns());
   const VectorView vector_bottom_part = view(vector, this->top_left_block.number_rows(), result.size());
   VectorView result_top_part = view(result, 0, this->top_left_block.number_rows());
   VectorView result_bottom_part = view(result, this->top_left_block.number_rows(), result.size());

   // carry out the matrix-vector products block-wise
   this->top_left_block.product(vector_top_part, result_top_part);
   this->top_right_block.product(vector_bottom_part, result_top_part);
   transpose(this->top_right_block).product(vector_top_part, result_bottom_part);
   this->bottom_right_block.product(vector_bottom_part, result_bottom_part);
}

#endif // UNO_SYMMETRIC2BY2BLOCKMATRIX_H