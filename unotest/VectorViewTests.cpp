// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "symbolic/VectorView.hpp"

using namespace uno;

TEST(VectorView, Size) {
   const std::vector<double> x{1., 2., 3., 100., 200., 300.};
   // vector view
   const size_t start_index = 2;
   const size_t end_index = 5;
   const auto x_view = view(x, start_index, end_index);
   ASSERT_EQ(x_view.size(), end_index - start_index);
}