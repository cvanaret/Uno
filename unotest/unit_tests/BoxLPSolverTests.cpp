// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "../DenseTestModel.hpp"
#include "ingredients/hessian_models/ZeroHessian.hpp"
#include "ingredients/inertia_correction_strategies/NoInertiaCorrection.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem_solvers/BoxLPSolver.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "tools/Statistics.hpp"

using namespace uno;

// Bug: BoxLPSolver only set the displacement of the multiplier of the active bound. The multiplier of the inactive bound
// kept its current value (its displacement was 0 instead of -current value), so that a variable that moves from its
// upper bound to its lower bound (or vice versa) ended up with nonzero multipliers on both bounds.
TEST(BoxLPSolver, MultipliersOfInactiveBoundsAreReset) {
   // min x0 - x1 + 0 x2  s.t.  0 <= x <= 1 (no general constraints)
   const DenseTestModel model("box_lp", 3, 0, ProblemType::LINEAR,
      [](const Vector<double>& x) { return x[0] - x[1]; },
      [](const Vector<double>&, Vector<double>& g) { g[0] = 1.; g[1] = -1.; g[2] = 0.; },
      [](const Vector<double>&, Vector<double>&) { },
      [](const Vector<double>&, double*) { },
      [](const Vector<double>&, double, const Vector<double>&, double*) { },
      {0., 0., 0.}, {1., 1., 1.}, {}, {}, {0.5, 0.5, 0.5});
   const OptimizationProblem problem(model, 100.);
   ZeroHessian hessian_model(model.number_variables);
   NoInertiaCorrection inertia_correction_strategy;
   const Subproblem subproblem(problem, hessian_model, inertia_correction_strategy);

   BoxLPSolver solver;
   solver.initialize_memory(subproblem);
   Iterate current_iterate(3, 0);
   current_iterate.primals = Vector<double>{0.5, 0.5, 0.5};
   // multipliers of a previous iterate, on the opposite bounds
   current_iterate.multipliers.lower_bounds = Vector<double>{0., 2., 4.};
   current_iterate.multipliers.upper_bounds = Vector<double>{-3., 0., -5.};
   Evaluations evaluations(model);
   Statistics statistics;
   const WarmstartInformation warmstart_information{};
   const Vector<double> initial_point(3);
   const Direction& direction = solver.solve(statistics, subproblem, current_iterate, 10., initial_point, evaluations,
      warmstart_information);
   ASSERT_EQ(direction.status, SubproblemStatus::OPTIMAL);

   // primal solution of the LP: lower bound for x0, upper bound for x1, x2 unchanged
   EXPECT_DOUBLE_EQ(direction.primals[0], -0.5);
   EXPECT_DOUBLE_EQ(direction.primals[1], 0.5);
   EXPECT_DOUBLE_EQ(direction.primals[2], 0.);

   // new multipliers = current multipliers + displacement: those of the LP solution
   const Vector<double> expected_lower_bound_multipliers{1., 0., 0.};
   const Vector<double> expected_upper_bound_multipliers{0., -1., 0.};
   for (size_t variable_index = 0; variable_index < 3; ++variable_index) {
      EXPECT_DOUBLE_EQ(current_iterate.multipliers.lower_bounds[variable_index] + direction.multipliers.lower_bounds[variable_index],
         expected_lower_bound_multipliers[variable_index]) << "variable " << variable_index;
      EXPECT_DOUBLE_EQ(current_iterate.multipliers.upper_bounds[variable_index] + direction.multipliers.upper_bounds[variable_index],
         expected_upper_bound_multipliers[variable_index]) << "variable " << variable_index;
   }
}