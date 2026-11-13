// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>
#include "linear_algebra/DenseMatrix.hpp"
#include "linear_algebra/Vector.hpp"
#include "linear_algebra/View.hpp"

using namespace uno;

TEST(View, Size) {
   const std::vector<double> x{1., 2., 3., 100., 200., 300.};
   // vector view
   const size_t start_index = 2;
   const size_t end_index = 5;
   const auto x_view = view(x, start_index, end_index);
   ASSERT_EQ(x_view.size(), end_index - start_index);
}
// Bug: y = x - z first copied x into y, then read z. With z = y, the original values of y were lost (result 0)
TEST(View, SubtractionWithAliasedRightOperand) {
   const Vector<double> x{1., 2.};
   Vector<double> y{3., 4.};
   y.view() = x.view() - y.view();
   EXPECT_EQ(y[0], -2.);
   EXPECT_EQ(y[1], -2.);
}

TEST(View, SubtractionWithBothOperandsAliased) {
   Vector<double> y{3., 4.};
   y.view() = y.view() - y.view();
   EXPECT_EQ(y[0], 0.);
   EXPECT_EQ(y[1], 0.);
}

TEST(View, LinearCombinationWithAliasedOperands) {
   const Vector<double> x{1., 2.};
   Vector<double> y{3., 4.};
   y.view() = x.view() + 2. * y.view(); // y = x + 2 y
   EXPECT_EQ(y[0], 7.);
   EXPECT_EQ(y[1], 10.);
   y.view() = x.view() - 2. * y.view(); // y = x - 2 y
   EXPECT_EQ(y[0], -13.);
   EXPECT_EQ(y[1], -18.);
   y.view() = y.view() + 0.5 * y.view(); // y = 1.5 y
   EXPECT_EQ(y[0], -19.5);
   EXPECT_EQ(y[1], -27.);
}

TEST(View, PartialOverlapIsRejected) {
   std::vector<double> data{1., 2., 3., 4.};
   auto first = view(data, 0, 3);
   auto second = view(data, 1, 4);
   EXPECT_THROW(first = second - first, std::invalid_argument);
   // the operand is left untouched
   EXPECT_EQ(data[0], 1.);
}

TEST(View, TransposedMatrixVectorProductInPlaceIsRejected) {
   DenseMatrix<double> A(2, 2);
   A.entry(0, 0) = 1.;
   A.entry(1, 1) = 1.;
   Vector<double> x{1., 2.};
   EXPECT_THROW(x.view() = transpose(A) * x.view(), std::invalid_argument);
}

// Bug: the triangular solves checked the solution against U but never checked the RHS, and passed the RHS size as the
// BLAS dimension: a larger RHS read and wrote out of bounds
TEST(View, TriangularSolveDimensionMismatchIsRejected) {
   DenseMatrix<double> U(2, 2);
   U.entry(0, 0) = 2.;
   U.entry(0, 1) = 1.;
   U.entry(1, 1) = 4.;
   const Vector<double> large_rhs{1., 2., 3.};
   Vector<double> solution(2);
   EXPECT_THROW(solution.view() = inverse(upper_triangular(U)) * large_rhs.view(), std::invalid_argument);
   EXPECT_THROW(solution.view() = transpose(inverse(upper_triangular(U))) * large_rhs.view(), std::invalid_argument);
   // a smaller RHS was read out of bounds by the copy into the solution
   const Vector<double> small_rhs{1.};
   EXPECT_THROW(solution.view() = inverse(upper_triangular(U)) * small_rhs.view(), std::invalid_argument);
   EXPECT_THROW(solution.view() = transpose(inverse(upper_triangular(U))) * small_rhs.view(), std::invalid_argument);

   // consistent dimensions: U x = (4, 8) gives x = (1, 2)
   const Vector<double> rhs{4., 8.};
   solution.view() = inverse(upper_triangular(U)) * rhs.view();
   EXPECT_DOUBLE_EQ(solution[0], 1.);
   EXPECT_DOUBLE_EQ(solution[1], 2.);
   // Uᵀ x = (2, 9) gives x = (1, 2)
   const Vector<double> rhs_transposed{2., 9.};
   solution.view() = transpose(inverse(upper_triangular(U))) * rhs_transposed.view();
   EXPECT_DOUBLE_EQ(solution[0], 1.);
   EXPECT_DOUBLE_EQ(solution[1], 2.);
}
