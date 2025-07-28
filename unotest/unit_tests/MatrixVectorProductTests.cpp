// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include <vector>
#include "linear_algebra/Vector.hpp"
#include "symbolic/MatrixVectorProduct.hpp"

using namespace uno;

TEST(MatrixVectorProduct, Test) {
   // (3, 7)
   // (7, 11)
   /*
   Vector<size_t> row_indices{0, 0, 1, 1};
   Vector<size_t> column_indices{0, 1, 0, 1};
   Vector<double> values{3., 7., 7., 11.};
   const std::vector<double> x{-2., 3.};
   const auto result = matrix * x;
   const std::vector<double> reference_result{15., 19.};
   for (size_t i: Range(result.size())) {
      ASSERT_EQ(result[i], reference_result[i]);
   }
   */
}
