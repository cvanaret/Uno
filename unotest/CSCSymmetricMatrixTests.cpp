// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project root for details.

#include <gtest/gtest.h>
#include "linear_algebra/CSCSymmetricMatrix.hpp"

TEST(CSCSymmetricMatrix, Identity) {
   CSCSymmetricMatrix m = CSCSymmetricMatrix::identity(3);
   // test that the entries are diagonal and equal to 1
   m.for_each([](int i, int j, double entry) {
      ASSERT_EQ(i, j);
      ASSERT_EQ(entry, 1.);
   });
}

// auxiliary function
CSCSymmetricMatrix generate_test_matrix(size_t padding_size) {
   const size_t dimension = 4;
   const size_t capacity = 4;
   CSCSymmetricMatrix m(dimension, capacity, padding_size);
   m.insert(10., 0, 0);
   m.finalize(0);
   m.insert(11., 0, 1);
   m.insert(12., 1, 1);
   m.finalize(1);
   m.finalize(2);
   m.insert(13., 1, 3);
   m.finalize(3);
   return m;
}

TEST(CSCSymmetricMatrix, Construction) {
   CSCSymmetricMatrix m = generate_test_matrix(1);

   ASSERT_EQ(m.column_start[0], 0);
   ASSERT_EQ(m.column_start[1], 1);
   ASSERT_EQ(m.column_start[2], 3);
   ASSERT_EQ(m.column_start[3], 3);
   ASSERT_EQ(m.column_start[4], 4);
}

TEST(CSCSymmetricMatrix, Regularization) {
   CSCSymmetricMatrix m = generate_test_matrix(1);
   // add 100*identity to m
   m.add_identity_multiple(100.);

   ASSERT_EQ(m.column_start[0], 0);
   ASSERT_EQ(m.column_start[1], 1);
   ASSERT_EQ(m.column_start[2], 3);
   ASSERT_EQ(m.column_start[3], 4);
   ASSERT_EQ(m.column_start[4], 6);
}

/*
TEST(CSCSymmetricMatrix, Filtering) {
   CSCSymmetricMatrix m = generate_test_matrix(0);
   std::cout << m;
   std::vector<int> active_set_variables{0, 2, 3};
   m.remove_variables(active_set_variables);
   std::cout << m;
}
*/