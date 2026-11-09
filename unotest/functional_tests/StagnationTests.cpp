// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Termination when the tolerances cannot be reached: once Uno has converged to the best possible accuracy, it accepts
// the same iterate with a zero step over and over. It must stop instead of iterating until the iteration limit.

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
   // GEMSEO's Power2 problem: min x0^2 + x1^2 + x2^2  s.t.  x0^3 >= 0.5, x1^3 >= 0.5, x2^3 = 0.9, -1 <= x <= 1
   DenseTestModel power2_model() {
      return DenseTestModel("power2", 3, 3, ProblemType::NONLINEAR,
         [](const Vector<double>& x) { return x[0]*x[0] + x[1]*x[1] + x[2]*x[2]; },
         [](const Vector<double>& x, Vector<double>& g) { g[0] = 2.*x[0]; g[1] = 2.*x[1]; g[2] = 2.*x[2]; },
         [](const Vector<double>& x, Vector<double>& c) { c[0] = x[0]*x[0]*x[0]; c[1] = x[1]*x[1]*x[1]; c[2] = x[2]*x[2]*x[2]; },
         [](const Vector<double>& x, double* J) {
            // dense Jacobian, row-major
            J[0] = 3.*x[0]*x[0]; J[1] = 0.; J[2] = 0.;
            J[3] = 0.; J[4] = 3.*x[1]*x[1]; J[5] = 0.;
            J[6] = 0.; J[7] = 0.; J[8] = 3.*x[2]*x[2];
         },
         [](const Vector<double>& x, double rho, const Vector<double>& y, double* H) {
            // dense lower triangle: (0,0), (1,0), (1,1), (2,0), (2,1), (2,2). Uno's convention: L = rho f - y^T c
            H[0] = 2.*rho - 6.*y[0]*x[0];
            H[1] = 0.; H[2] = 2.*rho - 6.*y[1]*x[1];
            H[3] = 0.; H[4] = 0.; H[5] = 2.*rho - 6.*y[2]*x[2];
         },
         {-1., -1., -1.}, {1., 1., 1.}, {0.5, 0.5, 0.9}, {Inf, Inf, 0.9}, {1., 1., 1.});
   }

   // tolerances that cannot be reached, acceptable termination disabled
   Options unreachable_tolerances_options(const std::string& preset) {
      Options options;
      DefaultOptions::load(options);
      Presets::set(options, preset);
      options.set_double("primal_tolerance", 0.);
      options.set_double("dual_tolerance", 0.);
      options.set_double("acceptable_tolerance", 0.);
      options.set_integer("acceptable_iteration_threshold", 0);
      options.set_integer("max_iterations", 1000);
      Logger::level = SILENT;
      return options;
   }

   void check_stagnation_termination(const Result& result, const std::string& label) {
      EXPECT_EQ(result.optimization_status, OptimizationStatus::SUCCESS) << label;
      EXPECT_EQ(result.solution_status, SolutionStatus::FEASIBLE_SMALL_STEP) << label;
      EXPECT_LT(result.number_iterations, 100) << label;
      const double x01 = std::cbrt(0.5), x2 = std::cbrt(0.9);
      EXPECT_NEAR(result.primal_solution[0], x01, 1e-6) << label;
      EXPECT_NEAR(result.primal_solution[1], x01, 1e-6) << label;
      EXPECT_NEAR(result.primal_solution[2], x2, 1e-6) << label;
      EXPECT_NEAR(result.solution_objective, 2.*x01*x01 + x2*x2, 1e-6) << label;
   }
} // namespace

// Bug: with zero tolerances and the acceptable test disabled, the IPOPT preset kept accepting the same iterate (up to
// round-off) until the iteration limit (or forever, in GEMSEO, whose limits do not count repeated evaluations)
TEST(Stagnation, InteriorPointStopsAtZeroStep) {
   if (SymmetricIndefiniteLinearSolverFactory::available_solvers().empty()) {
      GTEST_SKIP() << "the interior-point method requires a linear solver";
   }
   const DenseTestModel model = power2_model();
   Uno uno;
   const Result result = uno.solve(model, unreachable_tolerances_options("ipopt"));
   check_stagnation_termination(result, "ipopt");
}

// Bug: same with filter SQP, which alternated between two iterates that differ by round-off (f-type and h-type steps).
// With the line search, the tiny-step termination reported the converged point as INFEASIBLE_SMALL_STEP (a round-off
// infeasibility exceeds a zero primal tolerance)
TEST(Stagnation, FilterSQPStopsAtZeroStep) {
#ifndef HAS_BQPD
   GTEST_SKIP() << "the filtersqp preset requires BQPD";
#else
   for (const std::string globalization_mechanism: {"TR", "LS"}) {
      const DenseTestModel model = power2_model();
      Options options = unreachable_tolerances_options("filtersqp");
      options.set_string("globalization_mechanism", globalization_mechanism);
      Uno uno;
      const Result result = uno.solve(model, options);
      check_stagnation_termination(result, "filtersqp " + globalization_mechanism);
      EXPECT_LT(result.solution_stationarity, 1e-12) << globalization_mechanism;
      EXPECT_LT(result.solution_primal_feasibility, 1e-12) << globalization_mechanism;
   }
#endif
}

// with reachable tolerances, the termination is unchanged
TEST(Stagnation, ReachableTolerancesStillGiveKKTPoint) {
   if (SymmetricIndefiniteLinearSolverFactory::available_solvers().empty()) {
      GTEST_SKIP() << "the interior-point method requires a linear solver";
   }
   const DenseTestModel model = power2_model();
   Options options;
   DefaultOptions::load(options);
   Presets::set(options, "ipopt");
   Logger::level = SILENT;
   Uno uno;
   const Result result = uno.solve(model, options);
   EXPECT_EQ(result.optimization_status, OptimizationStatus::SUCCESS);
   EXPECT_EQ(result.solution_status, SolutionStatus::FEASIBLE_KKT_POINT);
}
