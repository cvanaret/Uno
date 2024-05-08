// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "tools/ChainCollection.hpp"
#include "tools/CollectionAdapter.hpp"
#include "tools/Range.hpp"
#include <vector>

TEST(ChainCollection, Size) {
   const std::vector<int> x{1, 2, 3};
   const std::vector<int> y{4, 5, 6};
   const auto chain = concatenate(CollectionAdapter(x), CollectionAdapter(y));
   /*chain.for_each([](size_t element) {
      std::cout << element << ' ';
   });
   std::cout << '\n';*/
   ASSERT_EQ(chain.size(), x.size() + y.size());
}

TEST(ChainCollection, Range) {
   const std::vector<size_t> x{5, 6, 7};
   const auto range = CollectionAdapter(Range(5));
   const auto chain = concatenate(range, CollectionAdapter(x));
   ASSERT_EQ(chain.size(), range.size() + x.size());
}

/*
TEST(ChainCollection, Iterators) {
   std::vector<size_t> x{5, 6, 7};
   const auto range = Range(100, 105);
   const auto chain = concatenate(CollectionAdapter(x), range);
   for (size_t index: chain) {
      std::cout << index << ' ';
   }
   std::cout << "\n";
}
*/