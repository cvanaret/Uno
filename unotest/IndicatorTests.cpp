// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include <vector>
#include "symbolic/Indicator.hpp"
#include "symbolic/Range.hpp"

TEST(Indicator, Cutoff) {
   const std::vector<double> x{1., 2., 3., 4., 5., 6., 7.};
   const auto indicator = (x < 4.5);
   const std::vector<double> reference_result{1., 1., 1., 1., 0., 0., 0.};
   for (size_t i: Range(x.size())) {
      ASSERT_EQ(indicator[i], reference_result[i]);
   }
}