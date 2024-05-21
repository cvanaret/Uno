// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "symbolic/CollectionAdapter.hpp"
#include <vector>

TEST(CollectionAdapter, Size) {
   const std::vector<int> x{1, 2, 3};
   const auto y = CollectionAdapter(x);
   ASSERT_EQ(y.size(), x.size());
}

TEST(CollectionAdapter, Iterator) {
   std::vector<size_t> x{5, 6, 7};
   const auto y = CollectionAdapter(x);
   for (const auto [index, value]: y) {
      ASSERT_EQ(value, x[index]);
   }
}