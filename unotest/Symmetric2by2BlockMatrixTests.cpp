// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "linear_algebra/COOSymmetricMatrix.hpp"
#include "linear_algebra/Symmetric2by2BlockMatrix.hpp"

COOSymmetricMatrix<double> identity(size_t dimension) {
   COOSymmetricMatrix<double> matrix(dimension, dimension, false);
   for (size_t index: Range(dimension)) {
      matrix.insert(1., index, index);
      matrix.finalize_column(index);
   }
   return matrix;
}

TEST(Symmetric2by2BlockMatrix, Correctness) {
   const size_t block_dimension = 2;
   const auto identity_block = identity(block_dimension);
   const auto zero_block = COOSymmetricMatrix<double>(block_dimension, block_dimension, false);

   const Symmetric2by2BlockMatrix block_matrix{identity_block, identity_block, identity_block};
   const std::vector<double> x{1., 2., 100., 200.};

   // matrix-vector product
   std::vector<double> result(x.size());
   initialize_vector(result, 0.);
   block_matrix.product(x, result);
   ASSERT_EQ(result[0], x[0] + x[2]);
   ASSERT_EQ(result[1], x[1] + x[3]);
}