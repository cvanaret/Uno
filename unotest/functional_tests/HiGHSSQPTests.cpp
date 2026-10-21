// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Regression tests of the inequality-constrained SQP path with the HiGHS QP solver

#include <gtest/gtest.h>
#include "../DenseTestModel.hpp"
#include "Uno.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "optimization/Result.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"
#include "tools/Logger.hpp"

using namespace uno;

namespace {
   // min -(x0^2 + x1^2) + 0.1 x0  s.t.  x0 + x1 <= 1.5,  0 <= x <= 1
   // The Lagrangian Hessian -2 I is indefinite: the QP is solvable only after a primal regularization of the Hessian.
   // Starting from (0.3, 0.4), SQP converges to the KKT point (0.5, 1) with f = -1.2.
   DenseTestModel nonconvex_model() {
      return DenseTestModel("nonconvex", 2, 1, ProblemType::NONLINEAR,
         [](const Vector<double>& x) { return -(x[0]*x[0] + x[1]*x[1]) + 0.1*x[0]; },
         [](const Vector<double>& x, Vector<double>& g) { g[0] = -2.*x[0] + 0.1; g[1] = -2.*x[1]; },
         [](const Vector<double>& x, Vector<double>& c) { c[0] = x[0] + x[1]; },
         [](const Vector<double>&, double* J) { J[0] = 1.; J[1] = 1.; },
         [](const Vector<double>&, double rho, const Vector<double>&, double* H) { H[0] = -2.*rho; H[1] = 0.; H[2] = -2.*rho; },
         {0., 0.}, {1., 1.}, {-Inf}, {1.5}, {0.3, 0.4}, {0});
   }
} // namespace

// Bug: HiGHSQuadraticProgram scattered the COO Hessian into HiGHS' CSC layout *before* calling the inertia correction,
// which assumes the COO layout (Hessian values followed by the regularization slots). The inertia was thus estimated
// on a scrambled matrix and the regularization values were written into arbitrary CSC slots: the Hessian passed to
// HiGHS remained indefinite and HiGHS aborted ("negative curvature").
TEST(HiGHSSQP, PrimalRegularizationOfIndefiniteHessian) {
   if (SymmetricIndefiniteLinearSolverFactory::available_solvers().empty()) {
      GTEST_SKIP() << "the primal inertia correction requires a linear solver";
   }
   const DenseTestModel model = nonconvex_model();
   Options options;
   DefaultOptions::load(options);
   Presets::set(options, "filtersqp");
   options.set_string("globalization_mechanism", "LS");
   options.set_string("inertia_correction_strategy", "primal");
   options.set_string("QP_solver", "HiGHS");
   Logger::level = SILENT;

   Uno uno;
   const Result result = uno.solve(model, options);
   ASSERT_EQ(result.optimization_status, OptimizationStatus::SUCCESS);
   ASSERT_EQ(result.solution_status, SolutionStatus::FEASIBLE_KKT_POINT);
   EXPECT_NEAR(result.primal_solution[0], 0.5, 1e-6);
   EXPECT_NEAR(result.primal_solution[1], 1., 1e-6);
   EXPECT_NEAR(result.solution_objective, -1.2, 1e-6);
}