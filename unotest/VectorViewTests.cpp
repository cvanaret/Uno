// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "linear_algebra/VectorView.hpp"

TEST(VectorView, size) {
   const std::vector<double> x{1., 2., 3., 100., 200., 300.};
   // vector view
   const size_t start_index = 2;
   const size_t end_index = 6;
   const VectorView x_view = VectorView(x, start_index, end_index);
   ASSERT_EQ(x_view.size(), end_index - start_index);
}