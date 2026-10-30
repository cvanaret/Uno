// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_REDUCTIONTESTPROBLEMS_H
#define UNO_REDUCTIONTESTPROBLEMS_H

#include "SparseTestModel.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   // min  (x0 - 3)^2 + (x1 - 1)^2 + (x2 - 3)^2 + x3^2 + x0 x1
   // s.t. c0: x0 + x1 <= 3         (linear, x1 fixed: bound x0 <= 1)
   //      c1: -2 x2 >= -4          (linear: bound x2 <= 2, negative coefficient)
   //      c2: x0 + x3 = 1          (linear, two free variables: kept)
   //      c3: x3^2 <= 4            (nonlinear with a single variable: kept)
   //      c4: 3 x2 >= -30          (linear: bound x2 >= -10, redundant with x2 >= 0)
   //      x1 = 2 (fixed), 0 <= x2 <= 10
   // solution x = (1, 2, 2, 0), f = 8, y = (-2, 1, 0, 0, 0), z_L = (0, 5, 0, 0), z_U = 0
   inline SparseTestProblem mixed_problem() {
      SparseTestProblem problem;
      problem.name = "mixed";
      problem.number_variables = 4;
      problem.number_constraints = 5;
      problem.objective = [](const Vector<double>& x) {
         return (x[0] - 3.)*(x[0] - 3.) + (x[1] - 1.)*(x[1] - 1.) + (x[2] - 3.)*(x[2] - 3.) + x[3]*x[3] + x[0]*x[1];
      };
      problem.gradient = [](const Vector<double>& x, Vector<double>& g) {
         g[0] = 2.*(x[0] - 3.) + x[1]; g[1] = 2.*(x[1] - 1.) + x[0]; g[2] = 2.*(x[2] - 3.); g[3] = 2.*x[3];
      };
      problem.constraints = [](const Vector<double>& x, Vector<double>& c) {
         c[0] = x[0] + x[1]; c[1] = -2.*x[2]; c[2] = x[0] + x[3]; c[3] = x[3]*x[3]; c[4] = 3.*x[2];
      };
      problem.jacobian_rows = {0, 0, 1, 2, 2, 3, 4};
      problem.jacobian_columns = {0, 1, 2, 0, 3, 3, 2};
      problem.jacobian = [](const Vector<double>& x, double* J) {
         J[0] = 1.; J[1] = 1.; J[2] = -2.; J[3] = 1.; J[4] = 1.; J[5] = 2.*x[3]; J[6] = 3.;
      };
      problem.hessian_rows = {0, 1, 1, 2, 3};
      problem.hessian_columns = {0, 0, 1, 2, 3};
      problem.hessian = [](const Vector<double>&, double rho, const Vector<double>& y, double* H) {
         H[0] = 2.*rho; H[1] = rho; H[2] = 2.*rho; H[3] = 2.*rho; H[4] = 2.*rho - 2.*y[3];
      };
      problem.variables_lower_bounds = {-Inf, 2., 0., -Inf};
      problem.variables_upper_bounds = {Inf, 2., 10., Inf};
      problem.constraints_lower_bounds = {-Inf, -4., 1., -Inf, -30.};
      problem.constraints_upper_bounds = {3., Inf, 1., 4., Inf};
      problem.initial_point = {0., 0., 1., 0.5}; // x1 is not at its fixed value
      problem.linear_constraints = {0, 1, 2, 4};
      return problem;
   }

   // elimination in rounds:
   // min  (x2 - 5)^2 + x1 x2 + x0 x1
   // s.t. c0: x0 + x1 = 3    (x0 fixed: x1 = 2 becomes fixed in round 1)
   //      c1: x1 + x2 <= 5   (x1 fixed in round 1: bound x2 <= 3)
   //      x0 = 1 (fixed)
   // solution x = (1, 2, 3), f = 12, y = (6, -2), z_L = 0, z_U = (-4, 0, 0)
   inline SparseTestProblem chained_problem() {
      SparseTestProblem problem;
      problem.name = "chained";
      problem.number_variables = 3;
      problem.number_constraints = 2;
      problem.objective = [](const Vector<double>& x) { return (x[2] - 5.)*(x[2] - 5.) + x[1]*x[2] + x[0]*x[1]; };
      problem.gradient = [](const Vector<double>& x, Vector<double>& g) {
         g[0] = x[1]; g[1] = x[2] + x[0]; g[2] = 2.*(x[2] - 5.) + x[1];
      };
      problem.constraints = [](const Vector<double>& x, Vector<double>& c) { c[0] = x[0] + x[1]; c[1] = x[1] + x[2]; };
      problem.jacobian_rows = {0, 0, 1, 1};
      problem.jacobian_columns = {0, 1, 1, 2};
      problem.jacobian = [](const Vector<double>&, double* J) { J[0] = 1.; J[1] = 1.; J[2] = 1.; J[3] = 1.; };
      problem.hessian_rows = {1, 2, 2};
      problem.hessian_columns = {0, 1, 2};
      problem.hessian = [](const Vector<double>&, double rho, const Vector<double>&, double* H) {
         H[0] = rho; H[1] = rho; H[2] = 2.*rho;
      };
      problem.variables_lower_bounds = {1., -Inf, -Inf};
      problem.variables_upper_bounds = {1., Inf, Inf};
      problem.constraints_lower_bounds = {3., -Inf};
      problem.constraints_upper_bounds = {3., 5.};
      problem.initial_point = {1., 0., 0.};
      problem.linear_constraints = {0, 1};
      return problem;
   }

   // all the constraints are declared nonlinear:
   // min  (x0 - 3)^2 + (x1 + 1)^2
   // s.t. c0: 2 x0 - x1 <= 10    (linear: detected)
   //      c1: x0^2 + x1 <= 5     (nonlinear)
   //      c2: x0 x2 <= 4         (linear in the reduced space: detected, then bound x0 <= 2)
   //      x2 = 2 (fixed)
   // solution x = (2, -1, 2), f = 1, y = (0, 0, -1), z_L = (0, 0, 2), z_U = 0
   inline SparseTestProblem undeclared_linear_problem() {
      SparseTestProblem problem;
      problem.name = "undeclared_linear";
      problem.number_variables = 3;
      problem.number_constraints = 3;
      problem.objective = [](const Vector<double>& x) { return (x[0] - 3.)*(x[0] - 3.) + (x[1] + 1.)*(x[1] + 1.); };
      problem.gradient = [](const Vector<double>& x, Vector<double>& g) { g[0] = 2.*(x[0] - 3.); g[1] = 2.*(x[1] + 1.); g[2] = 0.; };
      problem.constraints = [](const Vector<double>& x, Vector<double>& c) {
         c[0] = 2.*x[0] - x[1]; c[1] = x[0]*x[0] + x[1]; c[2] = x[0]*x[2];
      };
      problem.jacobian_rows = {0, 0, 1, 1, 2, 2};
      problem.jacobian_columns = {0, 1, 0, 1, 0, 2};
      problem.jacobian = [](const Vector<double>& x, double* J) {
         J[0] = 2.; J[1] = -1.; J[2] = 2.*x[0]; J[3] = 1.; J[4] = x[2]; J[5] = x[0];
      };
      problem.hessian_rows = {0, 1, 2};
      problem.hessian_columns = {0, 1, 0};
      problem.hessian = [](const Vector<double>&, double rho, const Vector<double>& y, double* H) {
         H[0] = 2.*rho - 2.*y[1]; H[1] = 2.*rho; H[2] = -y[2];
      };
      problem.variables_lower_bounds = {-Inf, -Inf, 2.};
      problem.variables_upper_bounds = {Inf, Inf, 2.};
      problem.constraints_lower_bounds = {-Inf, -Inf, -Inf};
      problem.constraints_upper_bounds = {10., 5., 4.};
      problem.initial_point = {0., 0., 2.};
      return problem;
   }
} // namespace

#endif // UNO_REDUCTIONTESTPROBLEMS_H
