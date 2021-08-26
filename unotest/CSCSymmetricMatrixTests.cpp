#include <gtest/gtest.h>
#include "linear_algebra/CSCSymmetricMatrix.hpp"

Level Logger::logger_level = INFO;

TEST(CSCSymmetricMatrix, Regularization) {
   const size_t dimension = 4;
   const size_t capacity = 4;
   const size_t padding_size = 1;
   CSCSymmetricMatrix m(dimension, capacity, padding_size);
   m.insert(10., 0, 0);
   m.finalize(0);
   m.insert(11., 0, 1);
   m.insert(12., 1, 1);
   m.finalize(1);
   m.finalize(2);
   m.insert(13., 1, 3);
   m.finalize(3);

   ASSERT_EQ(m.column_start[0], 0);
   ASSERT_EQ(m.column_start[1], 1);
   ASSERT_EQ(m.column_start[2], 3);
   ASSERT_EQ(m.column_start[3], 3);
   ASSERT_EQ(m.column_start[4], 4);

   // add 100*identity to m
   m.add_identity_multiple(100.);

   ASSERT_EQ(m.column_start[0], 0);
   ASSERT_EQ(m.column_start[1], 1);
   ASSERT_EQ(m.column_start[2], 3);
   ASSERT_EQ(m.column_start[3], 4);
   ASSERT_EQ(m.column_start[4], 6);
}

 void test_compressed_matrix() {

}
