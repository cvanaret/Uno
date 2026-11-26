// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Ingredients of the soft restoration phase of the line search (IPOPT, Wächter & Biegler 2006, Section 3.3):
// - the primal-dual error of the barrier problem (IpoptCalculatedQuantities::trial_primal_dual_system_error);
// - the acceptance test of the filter without predicted reduction (FilterLSAcceptor::CheckAcceptabilityOfTrialPoint(0)),
//   which must not modify the filter

#include <cmath>
#include <limits>
#include <gtest/gtest.h>
#include "ingredients/globalization_strategies/Decision.hpp"
#include "ingredients/globalization_strategies/ProgressMeasures.hpp"
#include "ingredients/globalization_strategies/switching_methods/filter_methods/WaechterFilterMethod.hpp"
#include "ingredients/inequality_handling_methods/interior_point_methods/InteriorPointParameters.hpp"
#include "ingredients/inequality_handling_methods/interior_point_methods/barrier_problems/PrimalDualInteriorPointProblem.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"
#include "tools/Statistics.hpp"
#include "BoundActiveModel.hpp"

using namespace uno;
using uno::test::BoundActiveModel;

// min (x0 + 1)^2 + (x1 - 2)^2 s.t. x0 + x1 >= 1, x0 >= 0, -1 <= x1 <= 1. Barrier problem in (x0, x1, s0), s0 >= 1
TEST(SoftRestoration, PrimalDualErrorOfBarrierProblem) {
   const BoundActiveModel model;
   const OptimizationProblem problem(model, 100.);
   const InteriorPointParameters parameters{0.99, 1e10, 0.25, 1e-2, 1e-2, 1e-5, 1.,
      std::pow(std::numeric_limits<double>::epsilon(), 0.75)};
   const double mu = 0.1;
   const PrimalDualInteriorPointProblem barrier_problem(problem, parameters, mu, 100.);
   ASSERT_EQ(barrier_problem.number_variables, 3u);
   ASSERT_EQ(barrier_problem.number_constraints, 1u);

   Iterate iterate(3, 1);
   iterate.primals[0] = 0.5;
   iterate.primals[1] = 0.5;
   iterate.primals[2] = 1.5;
   iterate.multipliers.constraints[0] = 0.2;
   iterate.multipliers.lower_bounds[0] = 1.;
   iterate.multipliers.lower_bounds[1] = 0.5;
   iterate.multipliers.lower_bounds[2] = 2.;
   iterate.multipliers.upper_bounds[1] = -0.25;
   Evaluations evaluations(model);
   evaluations.evaluate_objective_gradient(model, iterate.primals);
   evaluations.evaluate_jacobian(model, iterate.primals);

   // dual infeasibility: grad f - J^T y - z_L - z_U = (3 - 0.2 - 1, -3 - 0.2 - 0.5 + 0.25, 0.2 - 2), averaged over 3 entries
   const double dual_infeasibility = (1.8 + 3.45 + 1.8) / 3.;
   // primal infeasibility: c(x) - s = 1 - 1.5
   const double primal_infeasibility = 0.5;
   // mu-shifted complementarity: (1*0.5 - mu, 0.5*1.5 - mu, -0.25*(0.5 - 1) - mu, 2*(1.5 - 1) - mu), averaged over 4 bounds
   const double complementarity = (0.4 + 0.65 + 0.025 + 0.9) / 4.;
   EXPECT_NEAR(barrier_problem.primal_dual_system_error(iterate, evaluations),
      dual_infeasibility + primal_infeasibility + complementarity, 1e-14);
}

namespace {
   ProgressMeasures progress(double infeasibility, double objective) {
      return {infeasibility, [=](double) { return objective; }, 0., 1.};
   }

   Options ipopt_options() {
      Options options;
      DefaultOptions::load(options);
      Presets::set(options, "ipopt");
      return options;
   }
} // namespace

// the h-type acceptance test of the regular line search adds the current iterate to the filter; the acceptance test
// of the soft restoration phase leaves the filter untouched
TEST(SoftRestoration, AcceptanceTestDoesNotModifyTheFilter) {
   const Options options = ipopt_options();
   Statistics statistics{};
   statistics.add_column("Status", Statistics::string_width, 0);
   WaechterFilterMethod method(options);
   Iterate initial_iterate(1, 1);
   initial_iterate.progress = progress(1., 0.);
   method.initialize(statistics, initial_iterate);

   // halving the infeasibility is acceptable (h-type)
   EXPECT_TRUE(method.is_acceptable_without_predicted_reduction(progress(1., 0.), progress(0.5, 0.)));
   // (1, 0) is not in the filter: from (0.6, 1), the trial iterate (1, 0) reduces the merit and is acceptable
   EXPECT_TRUE(method.is_acceptable_without_predicted_reduction(progress(0.6, 1.), progress(1., 0.)));

   // the regular (h-type) acceptance test adds the current iterate (1, 0) to the filter
   EXPECT_TRUE(method.is_iterate_acceptable(statistics, progress(1., 0.), progress(0.5, 0.), progress(0., 0.), 1.).accepted);
   EXPECT_FALSE(method.is_acceptable_without_predicted_reduction(progress(0.6, 1.), progress(1., 0.)));
}

// as in the regular test, no h-type acceptance at a feasible current iterate
TEST(SoftRestoration, AcceptanceTestRejectsHTypeStepAtFeasibleIterate) {
   const Options options = ipopt_options();
   Statistics statistics{};
   statistics.add_column("Status", Statistics::string_width, 0);
   WaechterFilterMethod method(options);
   Iterate initial_iterate(1, 1);
   initial_iterate.progress = progress(1., 0.);
   method.initialize(statistics, initial_iterate);
   EXPECT_FALSE(method.is_acceptable_without_predicted_reduction(progress(0., 1.), progress(0., 0.5)));
}

// restoration phase entered from the soft restoration phase: the iterate at which the restoration phase started is not
// in the filter, and the exit test checks the sufficient reduction wrt this iterate (IPOPT, IsAcceptableToCurrentIterate)
TEST(SoftRestoration, AcceptabilityWrtRestorationReferenceIterate) {
   const Options options = ipopt_options();
   Statistics statistics{};
   WaechterFilterMethod method(options);
   Iterate initial_iterate(1, 1);
   initial_iterate.progress = progress(1., 0.);
   method.initialize(statistics, initial_iterate);
   const ProgressMeasures reference = progress(1., 0.);
   EXPECT_TRUE(method.is_acceptable_wrt_reference_iterate(reference, progress(0.5, 0.)));
   EXPECT_TRUE(method.is_acceptable_wrt_reference_iterate(reference, progress(1., -1.)));
   EXPECT_FALSE(method.is_acceptable_wrt_reference_iterate(reference, progress(1., 0.)));
}
