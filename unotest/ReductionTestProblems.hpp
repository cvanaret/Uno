// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_REDUCTIONTESTPROBLEMS_H
#define UNO_REDUCTIONTESTPROBLEMS_H

#include "SparseTestModel.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   // min  (x0 - 3)^2 + (x1 - 1)^2 + (x2 - 3)^2 + x3^2 + x0 x1
   // s.t. c0: x0 + x1 <= 3
   //      c1: -2 x2 >= -4
   //      c2: x0 + x3 = 1
   //      c3: x3^2 <= 4
   //      c4: x1^2 = 4        (fixed variable only, satisfied: removed)
   //      x1 = 2 (fixed), 0 <= x2 <= 10
   // solution x = (1, 2, 2, 0), f = 8, y = (-2, 1, 0, 0, 0), z_L = (0, 5, 0, 0), z_U = 0
   inline SparseTestProblem fixed_variable_problem() {
      SparseTestProblem problem;
      problem.name = "fixed_variable";
      problem.number_variables = 4;
      problem.number_constraints = 5;
      problem.objective = [](const Vector<double>& x) {
         return (x[0] - 3.)*(x[0] - 3.) + (x[1] - 1.)*(x[1] - 1.) + (x[2] - 3.)*(x[2] - 3.) + x[3]*x[3] + x[0]*x[1];
      };
      problem.gradient = [](const Vector<double>& x, Vector<double>& g) {
         g[0] = 2.*(x[0] - 3.) + x[1]; g[1] = 2.*(x[1] - 1.) + x[0]; g[2] = 2.*(x[2] - 3.); g[3] = 2.*x[3];
      };
      problem.constraints = [](const Vector<double>& x, Vector<double>& c) {
         c[0] = x[0] + x[1]; c[1] = -2.*x[2]; c[2] = x[0] + x[3]; c[3] = x[3]*x[3]; c[4] = x[1]*x[1];
      };
      problem.jacobian_rows = {0, 0, 1, 2, 2, 3, 4};
      problem.jacobian_columns = {0, 1, 2, 0, 3, 3, 1};
      problem.jacobian = [](const Vector<double>& x, double* J) {
         J[0] = 1.; J[1] = 1.; J[2] = -2.; J[3] = 1.; J[4] = 1.; J[5] = 2.*x[3]; J[6] = 2.*x[1];
      };
      problem.hessian_rows = {0, 1, 1, 2, 3};
      problem.hessian_columns = {0, 0, 1, 2, 3};
      problem.hessian = [](const Vector<double>&, double rho, const Vector<double>& y, double* H) {
         H[0] = 2.*rho; H[1] = rho; H[2] = 2.*rho - 2.*y[4]; H[3] = 2.*rho; H[4] = 2.*rho - 2.*y[3];
      };
      problem.variables_lower_bounds = {-Inf, 2., 0., -Inf};
      problem.variables_upper_bounds = {Inf, 2., 10., Inf};
      problem.constraints_lower_bounds = {-Inf, -4., 1., -Inf, 4.};
      problem.constraints_upper_bounds = {3., Inf, 1., 4., 4.};
      problem.initial_point = {0., 0., 1., 0.5}; // x1 is not at its fixed value
      problem.linear_constraints = {0, 1, 2};
      return problem;
   }

   // all the variables are fixed: x0 = 1, x1 = 2
   // min x0^2 + x1  s.t. c0: lower_bound <= x0 x1 <= 5 (constant: 2)
   // if lower_bound <= 2, c0 is removed: f = 3, y = 0, z_L = (2, 1), z_U = 0. Otherwise, c0 is kept and violated
   inline SparseTestProblem all_fixed_problem(double lower_bound) {
      SparseTestProblem problem;
      problem.name = "all_fixed";
      problem.number_variables = 2;
      problem.number_constraints = 1;
      problem.objective = [](const Vector<double>& x) { return x[0]*x[0] + x[1]; };
      problem.gradient = [](const Vector<double>& x, Vector<double>& g) { g[0] = 2.*x[0]; g[1] = 1.; };
      problem.constraints = [](const Vector<double>& x, Vector<double>& c) { c[0] = x[0]*x[1]; };
      problem.jacobian_rows = {0, 0};
      problem.jacobian_columns = {0, 1};
      problem.jacobian = [](const Vector<double>& x, double* J) { J[0] = x[1]; J[1] = x[0]; };
      problem.hessian_rows = {0, 1};
      problem.hessian_columns = {0, 0};
      problem.hessian = [](const Vector<double>&, double rho, const Vector<double>& y, double* H) { H[0] = 2.*rho; H[1] = -y[0]; };
      problem.variables_lower_bounds = {1., 2.};
      problem.variables_upper_bounds = {1., 2.};
      problem.constraints_lower_bounds = {lower_bound};
      problem.constraints_upper_bounds = {5.};
      problem.initial_point = {0., 0.};
      return problem;
   }
} // namespace

#endif // UNO_REDUCTIONTESTPROBLEMS_H
