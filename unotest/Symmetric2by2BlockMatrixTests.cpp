// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "linear_algebra/Symmetric2by2BlockMatrix.hpp"
#include "symbolic/MatrixExpression.hpp"

template <typename ElementType>
auto assemble_matrix(size_t block_dimension) {
   const Symmetric2by2BlockMatrix block_matrix(
         identity<ElementType>(block_dimension),       identity<ElementType>(block_dimension),
         /* identity<ElementType>(block_dimension) */  identity<ElementType>(block_dimension)
   );
   return block_matrix;
}

TEST(Symmetric2by2BlockMatrix, Correctness) {
   using T = double;

   const size_t block_dimension = 2;
   const Symmetric2by2BlockMatrix block_matrix = assemble_matrix<T>(block_dimension);
   const std::vector<T> x{1., 2., 100., 200.};

   // matrix-vector product
   std::vector<T> result(x.size(), 0.);
   block_matrix.product(x, result);
   const std::vector<T> reference_result{101., 202., 101., 202.};

   for (size_t index: Range(2*block_dimension)) {
      ASSERT_EQ(result[index], reference_result[index]);
   }
}

TEST(Symmetric2by2BlockMatrix, FloatCorrectness) {
   using T = float;

   const size_t block_dimension = 2;
   const Symmetric2by2BlockMatrix block_matrix = assemble_matrix<T>(block_dimension);
   const std::vector<T> x{1., 2., 100., 200.};

   // matrix-vector product
   std::vector<T> result(x.size(), 0.);
   block_matrix.product(x, result);
   const std::vector<T> reference_result{101., 202., 101., 202.};

   for (size_t index: Range(2*block_dimension)) {
      ASSERT_EQ(result[index], reference_result[index]);
   }
}