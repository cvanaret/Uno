// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "linear_algebra/COOSparseStorage.hpp"

using namespace uno;

const size_t n = 5;
const size_t nnz = 7;

COOSparseStorage<size_t, double> create_COO_storage() {
   COOSparseStorage<size_t, double> sparse_storage(n, nnz, false);
   sparse_storage.insert(2., 0, 0);
   sparse_storage.insert(3., 0, 1);
   sparse_storage.insert(4., 1, 2);
   sparse_storage.insert(6., 1, 4);
   sparse_storage.insert(1., 2, 2);
   sparse_storage.insert(5., 2, 3);
   sparse_storage.insert(1., 4, 4);
   return sparse_storage;
}

TEST(COOSparseStorage, NNZ) {
   const COOSparseStorage<size_t, double> sparse_storage = create_COO_storage();
   ASSERT_EQ(sparse_storage.number_nonzeros, nnz);
}