// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "linear_algebra/CSCSparseStorage.hpp"

using namespace uno;

const size_t n = 5;
const size_t nnz = 4;

CSCSparseStorage<size_t, double> create_CSC_storage() {
   CSCSparseStorage<size_t, double> sparse_storage(n, nnz, false);
   sparse_storage.insert(2., 0, 0);
   sparse_storage.finalize_column(0);
   sparse_storage.finalize_column(1);

   sparse_storage.insert(4., 1, 2);
   sparse_storage.insert(1., 2, 2);
   sparse_storage.finalize_column(2);

   sparse_storage.insert(5., 2, 3);
   sparse_storage.finalize_column(3);
   sparse_storage.finalize_column(4);
   return sparse_storage;
}

TEST(CSCSparseStorage, NNZ) {
   const CSCSparseStorage<size_t, double> sparse_storage = create_CSC_storage();
   ASSERT_EQ(sparse_storage.number_nonzeros, nnz);
}