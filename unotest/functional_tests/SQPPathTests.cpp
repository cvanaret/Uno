// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Regression tests of the SQP/SLP path (inequality-constrained subproblems, line search and trust region)

#include <cmath>
#include <string>
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
   Options sqp_options(const std::string& preset, const std::string& globalization_mechanism) {
      Options options;
      DefaultOptions::load(options);
      Presets::set(options, preset);
      options.set_string("globalization_mechanism", globalization_mechanism);
      options.set_string("QP_solver", "BQPD");
      options.set_string("LP_solver", "BQPD");
      Logger::level = SILENT;
      return options;
   }

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

// Bug: PrimalDualInertiaCorrection::regularize_hessian() threw "not implemented yet" (and its caller set up the
// augmented system instead of the Hessian system). Any inequality-constrained SQP run with
// inertia_correction_strategy = primal_dual stopped with an algorithmic error at the first iteration.
TEST(SQPPath, PrimalDualInertiaCorrectionOfQPHessian) {
   if (SymmetricIndefiniteLinearSolverFactory::available_solvers().empty()) {
      GTEST_SKIP() << "the inertia correction requires a linear solver";
   }
   for (const std::string globalization_mechanism: {"LS", "TR"}) {
      const DenseTestModel model = nonconvex_model();
      Options options = sqp_options("filtersqp", globalization_mechanism);
      options.set_string("inertia_correction_strategy", "primal_dual");

      Uno uno;
      const Result result = uno.solve(model, options);
      ASSERT_EQ(result.optimization_status, OptimizationStatus::SUCCESS) << globalization_mechanism;
      ASSERT_EQ(result.solution_status, SolutionStatus::FEASIBLE_KKT_POINT) << globalization_mechanism;
      EXPECT_NEAR(result.primal_solution[0], 0.5, 1e-6) << globalization_mechanism;
      EXPECT_NEAR(result.primal_solution[1], 1., 1e-6) << globalization_mechanism;
   }
}

// Bug: BacktrackingLineSearch only checked the INFEASIBLE and UNBOUNDED statuses of the optimality subproblem. A
// subproblem solver ERROR was treated as an optimal direction: the line search went on with whatever the solver left in
// the direction. Here (inconsistent variable bounds, a deterministic way to make BQPD fail) Uno accepted the bogus null
// step at every iteration until the iteration limit, and reported it as a regular run.
TEST(SQPPath, LineSearchHandlesSubproblemSolverError) {
   // min (x0 - 3)^2 + x1^2  s.t.  x0^2 + x1 <= 4,  -10 <= x0 <= 10,  1 <= x1 <= 0 (inconsistent)
   const DenseTestModel model("inconsistent_bounds", 2, 1, ProblemType::NONLINEAR,
      [](const Vector<double>& x) { return (x[0] - 3.)*(x[0] - 3.) + x[1]*x[1]; },
      [](const Vector<double>& x, Vector<double>& g) { g[0] = 2.*(x[0] - 3.); g[1] = 2.*x[1]; },
      [](const Vector<double>& x, Vector<double>& c) { c[0] = x[0]*x[0] + x[1]; },
      [](const Vector<double>& x, double* J) { J[0] = 2.*x[0]; J[1] = 1.; },
      [](const Vector<double>&, double rho, const Vector<double>& y, double* H) { H[0] = 2.*rho - 2.*y[0]; H[1] = 0.; H[2] = 2.*rho; },
      {-10., 1.}, {10., 0.}, {-Inf}, {4.}, {0.5, 0.5});
   Options options = sqp_options("filtersqp", "LS");
   options.set_integer("max_iterations", 50);

   Uno uno;
   const Result result = uno.solve(model, options);
   EXPECT_EQ(result.optimization_status, OptimizationStatus::ALGORITHMIC_ERROR);
   EXPECT_LE(result.number_iterations, 2);
}

// Bug: in the feasibility restoration phase, BacktrackingLineSearch did not check the status of the subproblem at all.
// When the linear constraints (not relaxed by the line-search restoration) are inconsistent, the infeasible
// restoration QP was used as a direction and accepted as a "0 primal step", with the multipliers of a failed QP solve.
TEST(SQPPath, LineSearchRejectsInfeasibleRestorationSubproblem) {
   // min x0^2 + x1^2  s.t.  x0 + x1 >= 3 (linear), x0 + x1 <= 1 (linear), x0^2 - x1 <= 0
   const DenseTestModel model("inconsistent_linear_constraints", 2, 3, ProblemType::NONLINEAR,
      [](const Vector<double>& x) { return x[0]*x[0] + x[1]*x[1]; },
      [](const Vector<double>& x, Vector<double>& g) { g[0] = 2.*x[0]; g[1] = 2.*x[1]; },
      [](const Vector<double>& x, Vector<double>& c) { c[0] = x[0] + x[1]; c[1] = x[0] + x[1]; c[2] = x[0]*x[0] - x[1]; },
      [](const Vector<double>& x, double* J) { J[0] = 1.; J[1] = 1.; J[2] = 1.; J[3] = 1.; J[4] = 2.*x[0]; J[5] = -1.; },
      [](const Vector<double>&, double rho, const Vector<double>& y, double* H) { H[0] = 2.*rho - 2.*y[2]; H[1] = 0.; H[2] = 2.*rho; },
      {-10., -10.}, {10., 10.}, {3., -Inf, -Inf}, {Inf, 1., 0.}, {0.5, 0.5}, {0, 1});
   const Options options = sqp_options("filtersqp", "LS");

   Uno uno;
   const Result result = uno.solve(model, options);
   EXPECT_EQ(result.optimization_status, OptimizationStatus::ALGORITHMIC_ERROR);
   // no iterate was accepted: the initial point is returned with its (zero) multipliers
   EXPECT_EQ(result.number_iterations, 1);
   for (size_t constraint_index = 0; constraint_index < 3; ++constraint_index) {
      EXPECT_EQ(result.constraint_dual_solution[constraint_index], 0.);
   }
}

// Bug: the small-step termination of the trust-region method was unreachable. The test "radius < TR_min_radius" was
// performed *before* the radius was decreased after a rejected step, and the decrease was immediately followed by
// "throw Small radius": Uno aborted with an algorithmic error even at a feasible point. Moreover, the termination test
// read the residuals of the *rejected* trial iterate, which are never computed (stale residuals of the Iterate object,
// i.e. those of the iterate preceding the current one).
// Uno should terminate with a feasible small step at x = 5. A large
// TR_min_radius triggers the small-radius logic deterministically at the second iteration:
// - iteration 1: the LP step to x = 10 is rejected (Armijo), the step to x = 5 is accepted
// - iteration 2: the step to x = 1 is rejected (no decrease) and the radius becomes 2 < TR_min_radius = 3
TEST(SQPPath, TrustRegionSmallStepTerminationUsesComputedResiduals) {
   // min (x - 3)^2  s.t.  x >= 1, from x = 0
   const DenseTestModel model("small_radius", 1, 1, ProblemType::NONLINEAR,
      [](const Vector<double>& x) { return (x[0] - 3.)*(x[0] - 3.); },
      [](const Vector<double>& x, Vector<double>& g) { g[0] = 2.*(x[0] - 3.); },
      [](const Vector<double>& x, Vector<double>& c) { c[0] = x[0]; },
      [](const Vector<double>&, double* J) { J[0] = 1.; },
      [](const Vector<double>&, double rho, const Vector<double>&, double* H) { H[0] = 2.*rho; },
      {-Inf}, {Inf}, {1.}, {Inf}, {0.});
   Options options = sqp_options("filterslp", "TR");
   options.set_double("TR_radius", 10.);
   options.set_double("TR_min_radius", 3.);

   Uno uno;
   const Result result = uno.solve(model, options);
   ASSERT_EQ(result.optimization_status, OptimizationStatus::SUCCESS);
   EXPECT_EQ(result.solution_status, SolutionStatus::FEASIBLE_SMALL_STEP);
   EXPECT_EQ(result.number_iterations, 2);
   // the reported point is the last accepted iterate, and its reported residuals are its own
   EXPECT_DOUBLE_EQ(result.primal_solution[0], 5.);
   EXPECT_DOUBLE_EQ(result.solution_objective, 4.);
   EXPECT_DOUBLE_EQ(result.solution_primal_feasibility, 0.);
}