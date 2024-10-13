// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "solvers/FortranCOOSparseStorage.hpp"

using namespace uno;

const size_t n = 2;
const size_t nnz = 3;

FortranCOOSparseStorage<size_t, double> create_Fortran_COO_storage() {
   FortranCOOSparseStorage<size_t, double> sparse_storage(n, nnz, false);
   sparse_storage.insert(4402., 0, 0);
   sparse_storage.insert(800., 0, 1);
   sparse_storage.insert(200., 1, 1);
   return sparse_storage;
}

// check that, in spite of their different internal representation, the indices are properly exposed
TEST(FortranCOOSparseStorage, NNZ) {
   const FortranCOOSparseStorage<size_t, double> sparse_storage = create_Fortran_COO_storage();
   size_t element_count = 0;
   for (const auto [row_index, column_index, element]: sparse_storage) {
      if (element_count == 0) {
         ASSERT_EQ(row_index, 0);
         ASSERT_EQ(column_index, 0);
      }
      if (element_count == 1) {
         ASSERT_EQ(row_index, 0);
         ASSERT_EQ(column_index, 1);
      }
      if (element_count == 2) {
         ASSERT_EQ(row_index, 1);
         ASSERT_EQ(column_index, 1);
      }
      element_count++;
   }
}