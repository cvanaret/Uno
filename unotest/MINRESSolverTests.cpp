// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "solvers/linear/iterative/MINRESSolver.hpp"
#include "linear_algebra/COOSymmetricMatrix.hpp"

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
   COOSymmetricMatrix<NumericalType> matrix(0, 0, false);
   std::vector<NumericalType> rhs{11., 18.};
   std::vector<NumericalType> result(dimension);
   MINRESSolver<NumericalType, decltype(linear_operator)> solver(linear_operator, dimension);

   solver.solve_indefinite_system(matrix, rhs, result);
   for (size_t i: Range(dimension)) {
      ASSERT_NEAR(result[i], reference_result[i], tolerance);
   }
}
