// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// residuals and local models of the primal-dual barrier problem (PrimalDualInteriorPointProblem)

#include <cmath>
#include <limits>
#include <gtest/gtest.h>
#include "BoundActiveModel.hpp"
#include "ingredients/globalization_strategies/PredictedReductionModels.hpp"
#include "ingredients/inequality_handling_methods/interior_point_methods/InteriorPointParameters.hpp"
#include "ingredients/inequality_handling_methods/interior_point_methods/barrier_problems/PrimalDualInteriorPointProblem.hpp"
#include "linear_algebra/Norm.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"

using namespace uno;
using uno::test::BoundActiveModel;

namespace {
   // PDIPM over min (x0 + 1)^2 + (x1 - 2)^2  s.t.  x0 + x1 >= 1,  x0 >= 0,  -1 <= x1 <= 1, with variables (x0, x1, s0)
   struct BarrierScenario {
      const BoundActiveModel model{};
      const OptimizationProblem problem{model, 100.};
      const InteriorPointParameters parameters{0.99, 1e10, 0.25, 1e-2, 1e-2, 1e-5, 1.,
         std::pow(std::numeric_limits<double>::epsilon(), 0.75)};
      double barrier_parameter{0.1};
      PrimalDualInteriorPointProblem barrier_problem{problem, parameters, barrier_parameter, 100.};

      [[nodiscard]] Iterate make_iterate(double x0, double x1, double s0) const {
         Iterate iterate(this->barrier_problem.number_variables, this->barrier_problem.number_constraints);
         iterate.primals[0] = x0;
         iterate.primals[1] = x1;
         iterate.primals[2] = s0;
         return iterate;
      }
   };
} // namespace

// c(x) = x0 + x1 = 1 satisfies the constraint c(x) >= 1, but the slack s0 = 3 lags behind: the point is not primal
// feasible for the slacked problem and its primal infeasibility must not vanish
TEST(BarrierProblem, PrimalInfeasibilityIncludesSlackResiduals) {
   BarrierScenario scenario;
   const Iterate iterate = scenario.make_iterate(0.5, 0.5, 3.);
   Evaluations evaluations(scenario.model);
   EXPECT_EQ(scenario.problem.primal_infeasibility(iterate, evaluations, Norm::INF), 0.);
   EXPECT_DOUBLE_EQ(scenario.barrier_problem.primal_infeasibility(iterate, evaluations, Norm::INF), 2.);
}

// consistent slack: the primal infeasibility is the constraint violation of the original problem
TEST(BarrierProblem, PrimalInfeasibilityWithConsistentSlack) {
   BarrierScenario scenario;
   const Iterate iterate = scenario.make_iterate(0.25, 0.25, 0.5); // c(x) = 0.5 < 1
   Evaluations evaluations(scenario.model);
   EXPECT_DOUBLE_EQ(scenario.barrier_problem.primal_infeasibility(iterate, evaluations, Norm::INF), 0.5);
}

// the predicted infeasibility reduction built at the current iterate must not change when the progress measures of a
// trial iterate are evaluated (the backtracking line search queries the model after evaluating trial iterates)
TEST(BarrierProblem, PredictedInfeasibilityReductionIsNotOverwritten) {
   BarrierScenario scenario;
   const Iterate current_iterate = scenario.make_iterate(0.5, 0.5, 3.);
   Evaluations current_evaluations(scenario.model);
   current_evaluations.evaluate_jacobian(scenario.model, current_iterate.primals);
   const Vector<double> direction{0.1, 0.2, -1.};
   const PredictedInfeasibilityReduction infeasibility_pred = scenario.barrier_problem.build_predicted_infeasibility_reduction(
      current_iterate, direction, Norm::L1, current_evaluations);
   const double reduction = infeasibility_pred(0.5);

   Iterate trial_iterate = scenario.make_iterate(10., -0.5, 1.5);
   Evaluations trial_evaluations(scenario.model);
   scenario.barrier_problem.set_infeasibility_measure(trial_iterate, trial_evaluations, Norm::L1);
   EXPECT_DOUBLE_EQ(infeasibility_pred(0.5), reduction);
   // c - s = -2 and J d = 0.1 + 0.2 + 1 (slack column -1): |-2 + 0.5 * 1.3| = 1.35
   EXPECT_DOUBLE_EQ(reduction, 2. - 1.35);
}
