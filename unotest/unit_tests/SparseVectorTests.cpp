// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "linear_algebra/SparseVector.hpp"

using namespace uno;

TEST(SparseVector, Empty) {
   SparseVector<double> x(2);
   // x is empty: its size (number of elements) is 0
   ASSERT_EQ(x.size(), 0);
}

TEST(SparseVector, Size) {
   const size_t capacity = 3;
   SparseVector<double> x(capacity);
   const double constant_term = 1.;
   x.insert(0, constant_term);
   x.insert(3, constant_term);
   x.insert(7, constant_term);
   ASSERT_EQ(x.size(), capacity);
}

TEST(SparseVector, AllIdentical) {
   const size_t capacity = 3;
   SparseVector<double> x(capacity);
   const double constant_term = 1.;
   x.insert(0, constant_term);
   x.insert(3, constant_term);
   x.insert(7, constant_term);
   for (const auto [index, entry]: x) {
      EXPECT_TRUE(index <= 7);
      ASSERT_EQ(entry, constant_term);
   }
}

TEST(SparseVector, Clear) {
   const size_t capacity = 2;
   SparseVector<double> x(capacity);
   x.insert(0, 1.);
   x.insert(3, 2.);
   // empty the vector
   x.clear();
   ASSERT_EQ(x.size(), 0);
}

TEST(SparseVector, InsertAfterClear) {
   const size_t capacity = 2;
   SparseVector<double> x(capacity);
   x.insert(0, 1.);
   x.insert(3, 2.);
   // empty the vector
   x.clear();
   // insert an element
   x.insert(7, 3.);
   ASSERT_EQ(x.size(), 1);
}
