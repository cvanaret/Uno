// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include <vector>
#include "linear_algebra/RectangularMatrix.hpp"
#include "symbolic/MatrixVectorProduct.hpp"

using namespace uno;

TEST(MatrixVectorProduct, Test) {
   // (3, 7)
   // (7, 11)
   const size_t dimension = 2;
   RectangularMatrix<double> matrix(dimension, dimension);
   matrix[0].insert(0, 3.);
   matrix[0].insert(1, 7.);
   matrix[1].insert(0, 7.);
   matrix[1].insert(1, 11.);
   const std::vector<double> x{-2., 3.};
   const auto result = matrix * x;
   const std::vector<double> reference_result{15., 19.};
   for (size_t i: Range(result.size())) {
      ASSERT_EQ(result[i], reference_result[i]);
   }
}
