// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "symbolic/Range.hpp"

using namespace uno;

TEST(Range, Size) {
   const size_t size = 5;
   const Range range(size);
   ASSERT_EQ(range.size(), size);
}

TEST(Range, Increasing) {
   const Range range(2, 5); // 2, 3, 4
   const std::vector<size_t> reference_result{2, 3, 4};
   size_t index = 0;
   for (size_t element: range) {
      ASSERT_EQ(element, reference_result[index]);
      index++;
   }
}

TEST(Range, Decreasing) {
   const BackwardRange range(5, 2); // 5, 4, 3
   const std::vector<size_t> reference_result{5, 4, 3};
   size_t index = 0;
   for (size_t element: range) {
      ASSERT_EQ(element, reference_result[index]);
      index++;
   }
}

TEST(Range, SumLoop) {
   const Range range(2, 5); // 2, 3, 4
   size_t sum = 0;
   for (size_t element: range) {
      sum += element;
   }
   ASSERT_EQ(sum, 9);
}

TEST(Range, BackwardSumLoop) {
   const BackwardRange range(5, 2); // 5, 4, 3
   size_t sum = 0;
   for (size_t element: range) {
      sum += element;
   }
   ASSERT_EQ(sum, 12);
}
