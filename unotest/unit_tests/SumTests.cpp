// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "linear_algebra/Vector.hpp"
#include "symbolic/Sum.hpp"
#include "symbolic/Range.hpp"

using namespace uno;

TEST(Sum, Test) {
   const Vector<double> x{100., 200., 300., 400.};
   const Vector<double> y{2., 3., 4., 5.};
   const Vector<double> reference_result{102., 203., 304., 405.};
   const auto sum = x + y;
   for (size_t i: Range(x.size())) {
      ASSERT_EQ(sum[i], reference_result[i]);
   }
}