// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "symbolic/Range.hpp"

TEST(Range, Size) {
   const size_t size = 5;
   const Range range(size);
   ASSERT_EQ(range.size(), size);
}

TEST(Range, SumLoop) {
   const Range range(2, 5);
   size_t sum = 0;
   for (size_t element: range) {
      sum += element;
   }
   ASSERT_EQ(sum, 9);
}

TEST(Range, BackwardSumLoop) {
   const BackwardRange range(5, 2);
   size_t sum = 0;
   for (size_t element: range) {
      sum += element;
   }
   ASSERT_EQ(sum, 12);
}

TEST(Range, SumRangeBasedLoop) {
   const Range range(2, 5);
   size_t sum = 0;
   for (const auto element: range) {
      sum += element;
   }
   ASSERT_EQ(sum, 9);
}

TEST(Range, BackwardSumForeach) {
   const BackwardRange range(5, 2);
   size_t sum = 0;
   for (const auto element: range) {
      sum += element;
   }
   ASSERT_EQ(sum, 12);
}