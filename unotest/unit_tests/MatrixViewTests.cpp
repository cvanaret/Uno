// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "linear_algebra/MatrixView.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"

using namespace uno;

TEST(MatrixView, SubmatrixDimension) {
   // create matrix
   const size_t dimension = 3;
   const size_t nnz = 6;
   SymmetricMatrix<size_t, double> matrix(dimension, nnz, false, "COO");

   // create 2x2 submatrix with 3 nonzeros
   const size_t submatrix_nnz = 3;
   auto submatrix = MatrixView<size_t, double>(matrix, 0, 0, 2, 2, submatrix_nnz);
   ASSERT_EQ(submatrix.get_number_rows(), 2);
   ASSERT_EQ(submatrix.get_number_columns(), 2);
}

TEST(MatrixView, WrongRow) {
   // create matrix
   const size_t dimension = 3;
   const size_t nnz = 6;
   SymmetricMatrix<size_t, double> matrix(dimension, nnz, false, "COO");

   // create 2x2 submatrix with 3 nonzeros
   const size_t submatrix_nnz = 3;
   auto submatrix = MatrixView<size_t, double>(matrix, 0, 0, 2, 2, submatrix_nnz);
   EXPECT_DEBUG_DEATH(submatrix.insert(3., 2, 0), "The row index exceeds the size of the matrix view");
}

TEST(MatrixView, WrongColumn) {
   // create matrix
   const size_t dimension = 3;
   const size_t nnz = 6;
   SymmetricMatrix<size_t, double> matrix(dimension, nnz, false, "COO");

   // create 2x2 submatrix with 3 nonzeros
   const size_t submatrix_nnz = 3;
   auto submatrix = MatrixView<size_t, double>(matrix, 0, 0, 2, 2, submatrix_nnz);
   EXPECT_DEBUG_DEATH(submatrix.insert(3., 0, 2), "The column index exceeds the size of the matrix view");
}