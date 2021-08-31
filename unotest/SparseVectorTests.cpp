#include <gtest/gtest.h>
#include "linear_algebra/SparseVector.hpp"

TEST(SparseVector, Empty) {
   SparseVector2<double> x(2);
   // x is empty: its size (number of elements) is 0
   ASSERT_EQ(x.size(), 0);
}

TEST(SparseVector, AllIdentical) {
   const size_t capacity = 3;
   SparseVector2<double> x(capacity);
   x.insert(0, 1.);
   x.insert(3, 1.);
   x.insert(7, 1.);
   ASSERT_EQ(x.size(), capacity);
}

//// auxiliary function
//CSCSymmetricMatrix generate_test_matrix(size_t padding_size) {
   //const size_t dimension = 4;
   //const size_t capacity = 4;
   //CSCSymmetricMatrix m(dimension, capacity, padding_size);
   //m.insert(10., 0, 0);
   //m.finalize(0);
   //m.insert(11., 0, 1);
   //m.insert(12., 1, 1);
   //m.finalize(1);
   //m.finalize(2);
   //m.insert(13., 1, 3);
   //m.finalize(3);
   //return m;
//}

//TEST(CSCSymmetricMatrix, Construction) {
   //CSCSymmetricMatrix m = generate_test_matrix(1);

   //ASSERT_EQ(m.column_start[0], 0);
   //ASSERT_EQ(m.column_start[1], 1);
   //ASSERT_EQ(m.column_start[2], 3);
   //ASSERT_EQ(m.column_start[3], 3);
   //ASSERT_EQ(m.column_start[4], 4);
//}

//TEST(CSCSymmetricMatrix, Regularization) {
   //CSCSymmetricMatrix m = generate_test_matrix(1);
   //// add 100*identity to m
   //m.add_identity_multiple(100.);

   //ASSERT_EQ(m.column_start[0], 0);
   //ASSERT_EQ(m.column_start[1], 1);
   //ASSERT_EQ(m.column_start[2], 3);
   //ASSERT_EQ(m.column_start[3], 4);
   //ASSERT_EQ(m.column_start[4], 6);
//}
