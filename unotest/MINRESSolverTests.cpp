// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "solvers/linear/iterative/MINRESSolver.hpp"
#include "linear_algebra/COOSymmetricMatrix.hpp"

using IndexType = size_t;
using NumericalType = double;

void linear_operator(const std::vector<NumericalType>& x, std::vector<NumericalType>& result) {
   result[0] = 1.*x[0] + 2.*x[1];
   result[1] = 2.*x[0] + 3.*x[1];
}

TEST(MINRESSolver, TwoDimensional) {
   const size_t dimension = 2;
   const std::vector<NumericalType> reference_result{3., 4.};
   const NumericalType tolerance{1e-6};
   // dummy matrix
   COOSymmetricMatrix<IndexType, NumericalType> matrix(0, 0, false);
   std::vector<NumericalType> rhs{11., 18.};
   std::vector<NumericalType> result(dimension);
   MINRESSolver<size_t, double, decltype(linear_operator)> solver(linear_operator, dimension);

   const bool from_scratch = true;
   solver.solve_indefinite_system(matrix, rhs, result, from_scratch);
   for (size_t i: Range(dimension)) {
      ASSERT_NEAR(result[i], reference_result[i], tolerance);
   }
}

void hs015_linear_operator(const std::vector<NumericalType>& x, std::vector<NumericalType>& result) {
   result[0] = 1.*x[0] + 1.*x[4] + 1.*x[5];
   result[1] = 1.*x[1] + -2.*x[4] + 2.*x[5];
   result[2] = 1.*x[2] + -1.*x[4];
   result[3] = 1.*x[3] + -1.*x[5];
   result[4] = 1.*x[0] + -2.*x[1] + -1.*x[2];
   result[5] = 1.*x[0] + 2.*x[1] + -1.*x[3];
}

TEST(MINRESSolver, Hs015LeastSquareDuals) {
   const size_t dimension = 6;
   // reference result of MA57
   const std::vector<NumericalType> reference_result{-33.6667, -2.77085, -28.125, -39.2084, -27.125, -38.2084};
   const NumericalType tolerance{1e-3};
   // dummy matrix
   COOSymmetricMatrix<IndexType, NumericalType> matrix(0, 0, false);
   std::vector<NumericalType> rhs{-99, -24.9377, -1, -1, 0, 0 };
   std::vector<NumericalType> result(dimension);
   MINRESSolver<size_t, double, decltype(hs015_linear_operator)> solver(hs015_linear_operator, dimension);

   const bool from_scratch = true;
   solver.solve_indefinite_system(matrix, rhs, result, from_scratch);
   for (size_t i: Range(dimension)) {
      ASSERT_NEAR(result[i], reference_result[i], tolerance);
   }
}