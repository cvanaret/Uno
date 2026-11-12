// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Functional tests of the L-BFGS Hessian model. The curvature safeguard of the (s, y) pairs is invariant to the scaling
// of the gradients: a scaled objective is solved as the original one

#include <cmath>
#include <stdexcept>
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
   // min σ [(x0 - 1)^2 + 10 (x1 - x0^2)^2]  s.t.  x0 + x1 <= 3
   DenseTestModel scaled_model(double sigma) {
      return DenseTestModel("scaled_rosenbrock", 2, 1, ProblemType::NONLINEAR,
         [=](const Vector<double>& x) { return sigma * (std::pow(x[0] - 1., 2) + 10. * std::pow(x[1] - x[0]*x[0], 2)); },
         [=](const Vector<double>& x, Vector<double>& g) {
            g[0] = sigma * (2.*(x[0] - 1.) - 40.*x[0]*(x[1] - x[0]*x[0]));
            g[1] = sigma * 20.*(x[1] - x[0]*x[0]);
         },
         [](const Vector<double>& x, Vector<double>& c) { c[0] = x[0] + x[1]; },
         [](const Vector<double>&, double* J) { J[0] = 1.; J[1] = 1.; },
         [](const Vector<double>&, double, const Vector<double>&, double*) { throw std::runtime_error("no Hessian"); },
         {-Inf, -Inf}, {Inf, Inf}, {-Inf}, {3.}, {-1., 1.});
   }
} // namespace

TEST(LBFGS, GradientScalingDoesNotAffectConvergence) {
   if (SymmetricIndefiniteLinearSolverFactory::available_solvers().empty()) {
      GTEST_SKIP() << "the interior-point method requires a linear solver";
   }
   for (const double sigma: {1., 1e-2, 1e4}) {
      const DenseTestModel model = scaled_model(sigma);
      Options options;
      DefaultOptions::load(options);
      Presets::set(options, "ipopt");
      options.set_string("hessian_model", "LBFGS");
      options.set_bool("use_function_scaling", false);
      options.set_integer("max_iterations", 500);
      Logger::level = SILENT;
      Uno uno;
      const Result result = uno.solve(model, options);
      ASSERT_EQ(result.optimization_status, OptimizationStatus::SUCCESS) << "sigma = " << sigma;
      // the unconstrained minimizer (1, 1) is feasible
      EXPECT_NEAR(result.primal_solution[0], 1., 1e-4) << "sigma = " << sigma;
      EXPECT_NEAR(result.primal_solution[1], 1., 1e-4) << "sigma = " << sigma;
   }
}

TEST(LSR1, NegativeCurvatureIsRegularizedBeforeComputingDirection) {
   if (SymmetricIndefiniteLinearSolverFactory::available_solvers().empty()) {
      GTEST_SKIP() << "requires a direct linear solver";
   }
   DenseTestModel model("concave_box", 1, 0, ProblemType::NONLINEAR,
      [](const Vector<double>& x) { return -x[0] * x[0] / 2.; },
      [](const Vector<double>& x, Vector<double>& g) { g[0] = -x[0]; },
      [](const Vector<double>&, Vector<double>&) { },
      [](const Vector<double>&, double*) { },
      [](const Vector<double>&, double w, const Vector<double>&, double* h) { h[0] = -w; },
      {-10.}, {10.}, {}, {}, {0.5});
   for (const std::string strategy: {"primal_dual", "primal"}) {
      Options options;
      DefaultOptions::load(options);
      Presets::set(options, "ipopt");
      options.set_string("hessian_model", "LSR1");
      options.set_string("inertia_correction_strategy", strategy);
      options.set_double("bound_relaxation_factor", 0.);
      options.set_integer("max_iterations", 100);
      Logger::level = SILENT;
      Uno solver;
      const Result result = solver.solve(model, options);
      ASSERT_EQ(result.optimization_status, OptimizationStatus::SUCCESS) << strategy;
      EXPECT_NEAR(result.primal_solution[0], 10., 1e-5) << strategy;
   }
}
