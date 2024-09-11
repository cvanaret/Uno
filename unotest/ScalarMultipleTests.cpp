// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include <vector>
#include "symbolic/ScalarMultiple.hpp"
#include "symbolic/Range.hpp"

using namespace uno;

TEST(ScalarMultiple, TimesThree) {
   const std::vector<double> x{1., 2., 3.};
   const std::vector<double> reference_result{3., 6., 9.};
   const auto y = 3.*x;
   for (size_t i: Range(y.size())) {
      ASSERT_EQ(y[i], reference_result[i]);
   }
}
