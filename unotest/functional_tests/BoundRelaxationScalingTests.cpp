// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Bug: the interior-point method relaxed the bounds of the *scaled* model (BoundRelaxedModel decorating ScaledModel).
// A relaxation of the scaled bound by delta allows a violation of the original constraint of up to delta / sigma, where
// sigma is the constraint scaling factor (as low as function_scaling_lower_bound = 1e-8). The iterates then converged
// to points that violate the original constraints by more than tight_unscaled_primal_infeasibility_tolerance, and the
// solver never terminated (CUTEst READING4). As in IPOPT, the original (unscaled) bounds must be relaxed.

#include <algorithm>
#include <vector>
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

TEST(BoundRelaxationScaling, OriginalConstraintsAreSatisfiedWithSmallConstraintScaling) {
   if (SymmetricIndefiniteLinearSolverFactory::available_solvers().empty()) {
      GTEST_SKIP() << "No linear solver available";
   }
   // min (x - 2)^2  s.t.  1e10 x <= 1e10: the solution x = 1 is on the bound of a constraint whose gradient is 1e10,
   // hence a constraint scaling factor of 100/1e10 = 1e-8
   constexpr double coefficient = 1e10;
   const DenseTestModel model("badly_scaled_constraint", 1, 1, ProblemType::NONLINEAR,
      [](const Vector<double>& x) { return (x[0] - 2.) * (x[0] - 2.); },
      [](const Vector<double>& x, Vector<double>& g) { g[0] = 2. * (x[0] - 2.); },
      [](const Vector<double>& x, Vector<double>& c) { c[0] = coefficient * x[0]; },
      [](const Vector<double>&, double* J) { J[0] = coefficient; },
      [](const Vector<double>&, double rho, const Vector<double>&, double* H) { H[0] = 2. * rho; },
      {-Inf}, {Inf}, {-Inf}, {coefficient}, {0.});
   Options options;
   DefaultOptions::load(options);
   Presets::set(options, "ipopt");
   options.set_integer("max_iterations", 200);
   Logger::level = SILENT;

   Uno uno;
   const Result result = uno.solve(model, options);
   ASSERT_EQ(result.optimization_status, OptimizationStatus::SUCCESS);
   ASSERT_EQ(result.solution_status, SolutionStatus::FEASIBLE_KKT_POINT);
   // violation of the original constraint
   const double violation = std::max(0., coefficient * result.primal_solution[0] - coefficient);
   EXPECT_LE(violation, options.get_double("tight_unscaled_primal_infeasibility_tolerance"));
   EXPECT_NEAR(result.primal_solution[0], 1., 1e-8);
}
