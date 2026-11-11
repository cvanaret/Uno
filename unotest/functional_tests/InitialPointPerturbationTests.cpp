// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Regression tests: when the derivatives cannot be evaluated at the initial point (e.g. CUTEst S365, whose gradient is
// not defined at the initial point), Uno retries from deterministic perturbations of the initial point instead of
// terminating with an evaluation error.

#include <gtest/gtest.h>
#include "../DenseTestModel.hpp"
#include "Uno.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "optimization/Result.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"
#include "tools/Logger.hpp"

using namespace uno;

namespace {
   // min (x0 - 1)^2 + (x1 - 2)^2  s.t.  x0 + x1 <= 1, whose gradient cannot be evaluated at the initial point (0, 0)
   class GradientFailsAtInitialPointModel: public DenseTestModel {
   public:
      GradientFailsAtInitialPointModel():
         DenseTestModel("gradient_fails_at_initial_point", 2, 1, ProblemType::NONLINEAR,
            [](const Vector<double>& x) { return (x[0] - 1.)*(x[0] - 1.) + (x[1] - 2.)*(x[1] - 2.); },
            [](const Vector<double>& x, Vector<double>& g) {
               if (x[0] == 0. && x[1] == 0.) {
                  throw GradientEvaluationError();
               }
               g[0] = 2.*(x[0] - 1.); g[1] = 2.*(x[1] - 2.);
            },
            [](const Vector<double>& x, Vector<double>& c) { c[0] = x[0] + x[1]; },
            [](const Vector<double>&, double* J) { J[0] = 1.; J[1] = 1.; },
            [](const Vector<double>&, double rho, const Vector<double>&, double* H) { H[0] = 2.*rho; H[1] = 0.; H[2] = 2.*rho; },
            {-Inf, -Inf}, {Inf, Inf}, {-Inf}, {1.}, {0., 0.}) { }
   };

   Options ipopt_options() {
      Options options;
      DefaultOptions::load(options);
      Presets::set(options, "ipopt");
      Logger::level = SILENT;
      return options;
   }
} // namespace

TEST(InitialPointPerturbation, SolvesFromPerturbedInitialPoint) {
   if (SymmetricIndefiniteLinearSolverFactory::available_solvers().empty()) {
      GTEST_SKIP() << "the interior-point method requires a linear solver";
   }
   const GradientFailsAtInitialPointModel model;
   Uno uno;
   const Result result = uno.solve(model, ipopt_options());
   ASSERT_EQ(result.optimization_status, OptimizationStatus::SUCCESS);
   EXPECT_NEAR(result.solution_objective, 2., 1e-6); // optimum (0, 1)
}

TEST(InitialPointPerturbation, DisabledPerturbationIsEvaluationError) {
   if (SymmetricIndefiniteLinearSolverFactory::available_solvers().empty()) {
      GTEST_SKIP() << "the interior-point method requires a linear solver";
   }
   const GradientFailsAtInitialPointModel model;
   Options options = ipopt_options();
   options.set_integer("initial_point_perturbation_attempts", 0);
   Uno uno;
   const Result result = uno.solve(model, options);
   EXPECT_EQ(result.optimization_status, OptimizationStatus::EVALUATION_ERROR);
}
