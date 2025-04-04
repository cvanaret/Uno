// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "linear_algebra/DenseMatrix.hpp"

using namespace uno;

TEST(DenseMatrix, Access) {
   DenseMatrix<double> M(5, 2);
   ASSERT_EQ(M.get(1, 1), 0.);
}