// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "linear_algebra/COOSymmetricMatrix.hpp"

const size_t n = 5;
const size_t nnz = 7;

COOSymmetricMatrix<double> create_COO_matrix() {
   COOSymmetricMatrix<double> matrix(n, nnz, false);
   matrix.insert(2., 0, 0);
   matrix.insert(3., 0, 1);
   matrix.insert(4., 1, 2);
   matrix.insert(6., 1, 4);
   matrix.insert(1., 2, 2);
   matrix.insert(5., 2, 3);
   matrix.insert(1., 4, 4);
   return matrix;
}

TEST(COOSymmetricMatrix, NNZ) {
   const COOSymmetricMatrix<double> matrix = create_COO_matrix();
   ASSERT_EQ(matrix.number_nonzeros, nnz);
}

/*
TEST(COOSymmetricMatrix, Iterator) {
   const COOSymmetricMatrix<double> matrix = create_COO_matrix();
   for (const auto [i, j, Mij]: matrix) {
      std::cout << "COO(" << i << ", " << j << ") = " << Mij << '\n';
   }
}
*/