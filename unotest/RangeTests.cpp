// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "tools/Range.hpp"

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
   const Range<BACKWARD> range(5, 2);
   size_t sum = 0;
   for (size_t element: range) {
      sum += element;
   }
   ASSERT_EQ(sum, 12);
}

TEST(Range, SumForeach) {
   const Range range(2, 5);
   size_t sum = 0;
   range.for_each([&](size_t, size_t element) {
      sum += element;
   });
   ASSERT_EQ(sum, 9);
}

TEST(Range, BackwardSumForeach) {
   const Range<BACKWARD> range(5, 2);
   size_t sum = 0;
   range.for_each([&](size_t, size_t element) {
      sum += element;
   });
   ASSERT_EQ(sum, 12);
}