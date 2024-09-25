// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "linear_algebra/Vector.hpp"

using namespace uno;

TEST(Vector, Empty) {
   Vector<double> x(0);
   // x is empty: its size (number of elements) is 0
   ASSERT_EQ(x.size(), 0);
}

TEST(Vector, Size) {
   const size_t capacity = 3;
   Vector<double> x{1., 1., 1.};
   ASSERT_EQ(x.size(), capacity);
}

TEST(Vector, AllIdentical) {
   const double constant_term = 1.;
   Vector<double> x{constant_term, constant_term, constant_term};
   for (const auto element: x) {
      ASSERT_EQ(element, constant_term);
   }
}

TEST(Vector, MoveConstructor) {
   const double constant_term = 1.;
   Vector<double> x{constant_term, constant_term, constant_term};
   Vector<double> y = std::move(x);

   for (const auto element: y) {
      ASSERT_EQ(element, constant_term);
   }
}

TEST(Vector, MoveAssignment) {
   const double constant_term = 1.;
   Vector<double> x{constant_term, constant_term, constant_term};
   Vector<double> y{2.*constant_term, 2.*constant_term, 2.*constant_term};
   y = std::move(x);

   for (const auto element: y) {
      ASSERT_EQ(element, constant_term);
   }
}

TEST(Vector, Swap) {
   const double constant_term = 1.;
   Vector<double> x{constant_term, constant_term, constant_term};
   Vector<double> y{2.*constant_term, 2.*constant_term, 2.*constant_term};
   std::swap(x, y);

   for (const auto element: x) {
      ASSERT_EQ(element, 2.*constant_term);
   }
   for (const auto element: y) {
      ASSERT_EQ(element, constant_term);
   }
}
