// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "linear_algebra/CSCSymmetricMatrix.hpp"

const size_t n = 5;
const size_t nnz = 4;

CSCSymmetricMatrix<size_t, double> create_CSC_matrix() {
   CSCSymmetricMatrix<size_t, double> matrix(n, nnz, false);
   matrix.insert(2., 0, 0);
   matrix.finalize_column(0);
   matrix.finalize_column(1);

   matrix.insert(4., 1, 2);
   matrix.insert(1., 2, 2);
   matrix.finalize_column(2);

   matrix.insert(5., 2, 3);
   matrix.finalize_column(3);
   matrix.finalize_column(4);
   return matrix;
}

TEST(CSCSymmetricMatrix, NNZ) {
   const CSCSymmetricMatrix<size_t, double> matrix = create_CSC_matrix();
   ASSERT_EQ(matrix.number_nonzeros, nnz);
}

/*
TEST(CSCSymmetricMatrix, Iterator) {
   const CSCSymmetricMatrix<double> matrix = create_CSC_matrix();
   for (const auto [i, j, Mij]: matrix) {
      std::cout << "CSC(" << i << ", " << j << ") = " << Mij << '\n';
   }
}
*/