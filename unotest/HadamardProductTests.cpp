// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include <vector>
#include "symbolic/Expression.hpp"

TEST(HadamardProduct, Masking) {
   const std::vector<double> mask{0., 1., 1., 0., 1.};
   const std::vector<double> x{100., 200., 300., 400., 500.};
   const auto hadamard_product = hadamard(mask, x);
   const std::vector<double> reference_result{0., 200., 300., 0., 500.};
   for (size_t i: Range(x.size())) {
      ASSERT_EQ(hadamard_product[i], reference_result[i]);
   }
}

TEST(HadamardProduct, MaskingCombination) {
   const std::vector<double> mask1{0., 0., 1., 0., 1.};
   const std::vector<double> x{100., 200., 300., 400., 500.};

   const std::vector<double> mask2{1., 0., 1., 1., 0.};
   const std::vector<double> y{1000., 2000., 3000., 4000., 5000.};
   const auto hadamard_product = hadamard(mask1, x) + hadamard(mask2, y);
   const std::vector<double> reference_result{1000., 0., 3300., 4000., 500.};
   for (size_t i: Range(x.size())) {
      ASSERT_EQ(hadamard_product[i], reference_result[i]);
   }
}