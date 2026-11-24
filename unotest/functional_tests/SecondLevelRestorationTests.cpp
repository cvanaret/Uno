// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Second-level restoration (IPOPT, RestoRestorationPhase): when the backtracking line search fails in the feasibility
// restoration phase, the primals are kept and the elastic variables are reset to their closed-form values, which
// reduces the infeasibility of the restoration problem. Without it, the line search failure in the restoration phase
// was an algorithmic error

#include <cmath>
#include <gtest/gtest.h>
#include "Uno.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "optimization/Result.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"
#include "tools/Logger.hpp"
#include "../DenseTestModel.hpp"

using namespace uno;

namespace {
   // min x0 + x1 s.t. x0^2 + x1^2 + 1 = 0, x0^3 - a x1 = 0: the first constraint can't be satisfied, the expected outcome
   // is an infeasible stationary point (the minimizer of the constraint violation, x = 0). From these starting points, the
   // line search fails once in the restoration phase (after about 200 iterations)
   Result solve_infeasible_problem(double a, double x0, double x1) {
      const DenseTestModel model("infeasible_circle", 2, 2, ProblemType::NONLINEAR,
         [](const Vector<double>& x) { return x[0] + x[1]; },
         [](const Vector<double>&, Vector<double>& gradient) { gradient[0] = 1.; gradient[1] = 1.; },
         [=](const Vector<double>& x, Vector<double>& constraints) {
            constraints[0] = x[0] * x[0] + x[1] * x[1] + 1.;
            constraints[1] = x[0] * x[0] * x[0] - a * x[1];
         },
         [=](const Vector<double>& x, double* jacobian) {
            jacobian[0] = 2. * x[0]; jacobian[1] = 2. * x[1];
            jacobian[2] = 3. * x[0] * x[0]; jacobian[3] = -a;
         },
         // Lagrangian rho f - y^T c
         [](const Vector<double>& x, double /*objective_multiplier*/, const Vector<double>& multipliers, double* hessian) {
            hessian[0] = -2. * multipliers[0] - 6. * multipliers[1] * x[0];
            hessian[1] = 0.;
            hessian[2] = -2. * multipliers[0];
         },
         {-INFINITY, -INFINITY}, {INFINITY, INFINITY}, {0., 0.}, {0., 0.}, {x0, x1});
      Options options;
      DefaultOptions::load(options);
      Presets::set(options, "ipopt");
      options.set_integer("max_iterations", 500);
      Logger::level = SILENT;
      Uno uno;
      return uno.solve(model, options);
   }
} // namespace

TEST(SecondLevelRestoration, LineSearchFailureInRestorationPhaseIsRecovered) {
   if (SymmetricIndefiniteLinearSolverFactory::available_solvers().empty()) {
      GTEST_SKIP() << "no linear solver available";
   }
   const Result result = solve_infeasible_problem(1., 3., -2.);
   EXPECT_EQ(result.optimization_status, OptimizationStatus::SUCCESS);
   EXPECT_EQ(result.solution_status, SolutionStatus::INFEASIBLE_STATIONARY_POINT);
}

TEST(SecondLevelRestoration, LineSearchFailureInRestorationPhaseIsRecoveredOtherStartingPoint) {
   if (SymmetricIndefiniteLinearSolverFactory::available_solvers().empty()) {
      GTEST_SKIP() << "no linear solver available";
   }
   const Result result = solve_infeasible_problem(3., 1.5, -0.5);
   EXPECT_EQ(result.optimization_status, OptimizationStatus::SUCCESS);
   EXPECT_EQ(result.solution_status, SolutionStatus::INFEASIBLE_STATIONARY_POINT);
}
