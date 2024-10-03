// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include <vector>
#include "symbolic/Expression.hpp"

using namespace uno;

TEST(Sum, Test) {
   const std::vector<double> x{100., 200., 300., 400.};
   const std::vector<double> y{2., 3., 4., 5.};
   const std::vector<double> reference_result{102., 203., 304., 405.};
   const auto sum = x + y;
   for (size_t i: Range(x.size())) {
      ASSERT_EQ(sum[i], reference_result[i]);
   }
}