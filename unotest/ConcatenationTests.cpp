// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "symbolic/Concatenation.hpp"
#include "symbolic/CollectionAdapter.hpp"
#include "symbolic/Range.hpp"
#include <vector>

using namespace uno;

TEST(Concatenation, Size) {
   const std::vector<int> x{1, 2, 3};
   const std::vector<int> y{4, 5, 6};
   const auto chain = concatenate(CollectionAdapter(x), CollectionAdapter(y));
   /*chain.for_each([](size_t element) {
      std::cout << element << ' ';
   });
   std::cout << '\n';*/
   ASSERT_EQ(chain.size(), x.size() + y.size());
}

TEST(Concatenation, Range) {
   const std::vector<size_t> x{5, 6, 7};
   const auto range = Range(5);
   const auto chain = concatenate(range, CollectionAdapter(x));
   ASSERT_EQ(chain.size(), range.size() + x.size());
}

TEST(Concatenation, Iterator) {
   const std::vector<size_t> x{5, 6, 7};
   const auto range = Range(100, 105);
   const auto chain = concatenate(CollectionAdapter(x), range);
   const std::vector<size_t> reference_result{5, 6, 7, 100, 101, 102, 103, 104};
   size_t index = 0;
   for (const size_t value: chain) {
      ASSERT_EQ(value, reference_result[index]);
      index++;
   }
   std::cout << "\n";
}
