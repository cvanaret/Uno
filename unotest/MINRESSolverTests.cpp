// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "solvers/linear/iterative/MINRESSolver.hpp"
#include "linear_algebra/COOSymmetricMatrix.hpp"

TEST(MINRESSolver, TwoDimensional) {
   using NumericalType = double;
   const size_t dimension = 2;
   constexpr auto linear_operator = [=](const std::vector<NumericalType>& x, std::vector<NumericalType>& result) {
      result[0] = 1.*x[0] + 2.*x[1];
      result[1] = 2.*x[0] + 3.*x[1];
   };

   COOSymmetricMatrix<NumericalType> matrix(0, 0, false);
   std::vector<NumericalType> rhs(dimension);
   rhs[0] = 11.;
   rhs[1] = 18.;
   std::vector<NumericalType> result(dimension);
   MINRESSolver<NumericalType, decltype(linear_operator)> solver(linear_operator, dimension);

   solver.solve_indefinite_system(matrix, rhs, result);
}
