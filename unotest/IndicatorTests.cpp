// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include <vector>
#include "symbolic/Indicator.hpp"
#include "linear_algebra/Vector.hpp"

TEST(Indicator, Test) {
   const std::vector<double> x{1., 2., 3., 4., 5., 6., 7.};
   const auto indicator = (x < 4.5);
   //print_vector(std::cout, indicator);
}