// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "linear_algebra/DenseMatrix.hpp"
#include "linear_algebra/Vector.hpp"

using namespace uno;

TEST(DenseMatrix, Access) {
   DenseMatrix<double> matrix(5, 2);
   ASSERT_EQ(matrix.entry(1, 1), 0.);
}

DenseMatrix<double> get_matrix() {
   const size_t dimension = 3;
   // lower triangular matrix expressed in column-major order
   DenseMatrix<double> matrix(dimension, dimension);
   matrix.entry(0, 0) = 4.;
   matrix.entry(1, 0) = 12.;
   matrix.entry(2, 0) = -16.;
   matrix.entry(1, 1) = 37.;
   matrix.entry(2, 1) = -43.;
   matrix.entry(2, 2) = 98.;
   /*
      4., 12., -16.,
      0., 37., -43.,
      0., 0., 98.
   */
   return matrix;
}

TEST(DenseMatrix, Access2) {
   const DenseMatrix<double> matrix = get_matrix();
   ASSERT_EQ(matrix.entry(2, 1), -43.);
}

TEST(DenseMatrix, ColumnSize) {
   const DenseMatrix<double> matrix = get_matrix();
   const auto second_column = matrix.column(1);
   ASSERT_EQ(second_column.size(), 3);
}

TEST(DenseMatrix, ColumnAccess) {
   const DenseMatrix<double> matrix = get_matrix();
   const auto second_column = matrix.column(1);
   ASSERT_EQ(second_column[1], 37.);
}

TEST(DenseMatrix, ColumnMagnitude) {
   const DenseMatrix<double> matrix = get_matrix();
   const auto second_column = matrix.column(1);
   const double squared_magnitude = dot(second_column, second_column);
   ASSERT_EQ(squared_magnitude, 3218.);
}

/*
TEST(DenseMatrix, ColumnOverwrite) {
   const DenseMatrix<double> matrix = get_matrix();
   auto second_column = matrix.column(1);
   const Vector<double> x{1., 2., 3.};
   // overwrite the second column with x
   second_column = x;
   ASSERT_EQ(matrix.entry(0, 1), 1.);
   ASSERT_EQ(matrix.entry(1, 1), 2.);
   ASSERT_EQ(matrix.entry(2, 1), 3.);
}
*/