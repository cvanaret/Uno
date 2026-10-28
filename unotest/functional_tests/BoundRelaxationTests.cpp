// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// IPOPT-style dynamic bound relaxation in PrimalDualInteriorPointProblem: when the iterate gets too close to (or
// crosses) a bound, the bound is relaxed and the iterate is left untouched

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
#include <gtest/gtest.h>
#include "Uno.hpp"
#include "ingredients/inequality_handling_methods/interior_point_methods/InteriorPointParameters.hpp"
#include "ingredients/inequality_handling_methods/interior_point_methods/barrier_problems/PrimalDualInteriorPointProblem.hpp"
#include "model/Model.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "optimization/Result.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"
#include "symbolic/IntegerRange.hpp"
#include "symbolic/Range.hpp"
#include "tools/Logger.hpp"
#include "tools/UserCallbacks.hpp"
#include "BoundActiveModel.hpp"

using namespace uno;
using uno::test::BoundActiveModel;

namespace {
   const double slack_move = std::pow(std::numeric_limits<double>::epsilon(), 0.75);

   InteriorPointParameters make_parameters() {
      return {0.99, 1e10, 0.25, 1e-2, 1e-2, 1e-5, 1., slack_move};
   }

   // optimality-like scenario: PDIPM over the model, with variables (x0, x1, s0)
   struct BarrierScenario {
      explicit BarrierScenario(double mu): barrier_parameter(mu) { }

      const BoundActiveModel model{};
      const OptimizationProblem problem{model, 100.};
      const InteriorPointParameters parameters{make_parameters()};
      double barrier_parameter;
      PrimalDualInteriorPointProblem barrier_problem{problem, parameters, barrier_parameter, 100.};

      Iterate make_iterate(double x0, double x1, double s0, double z_lower = 1., double z_upper = -1.) const {
         Iterate iterate(this->barrier_problem.number_variables, this->barrier_problem.number_constraints);
         iterate.primals[0] = x0;
         iterate.primals[1] = x1;
         iterate.primals[2] = s0;
         for (size_t variable_index: Range(this->barrier_problem.number_variables)) {
            iterate.multipliers.lower_bounds[variable_index] = z_lower;
            iterate.multipliers.upper_bounds[variable_index] = (variable_index == 1) ? z_upper : 0.;
         }
         return iterate;
      }

      [[nodiscard]] double barrier_measure(Iterate& iterate) const {
         this->barrier_problem.set_auxiliary_measure(iterate);
         return iterate.progress.auxiliary;
      }
   };
} // namespace

TEST(BoundRelaxation, NoRelaxationAwayFromBounds) {
   BarrierScenario scenario(0.1);
   Iterate iterate = scenario.make_iterate(0.5, 0.5, 1.5);
   Evaluations evaluations(scenario.model);
   scenario.barrier_problem.postprocess_iterate(iterate, evaluations);
   EXPECT_EQ(scenario.barrier_problem.get_number_bound_relaxations(), 0u);
   EXPECT_EQ(scenario.barrier_problem.get_internal_variables_lower_bounds()[0], 0.);
   EXPECT_EQ(scenario.barrier_problem.get_internal_variables_upper_bounds()[1], 1.);
   EXPECT_EQ(scenario.barrier_problem.get_internal_variables_lower_bounds()[2], 1.);
}

TEST(BoundRelaxation, IterateOnLowerBoundRelaxesBoundNotIterate) {
   BarrierScenario scenario(0.1);
   Iterate iterate = scenario.make_iterate(0., 0.5, 1.5); // x0 exactly on its bound
   Evaluations evaluations(scenario.model);
   EXPECT_TRUE(std::isinf(scenario.barrier_measure(iterate)));

   scenario.barrier_problem.postprocess_iterate(iterate, evaluations);
   EXPECT_EQ(iterate.primals[0], 0.); // the iterate is not moved
   EXPECT_EQ(scenario.barrier_problem.get_number_bound_relaxations(), 1u);
   // mu/z = 0.1 > cap: the new slack is slack_move*max(1, |0|)
   EXPECT_DOUBLE_EQ(scenario.barrier_problem.get_internal_variables_lower_bounds()[0], -slack_move);
   EXPECT_TRUE(std::isfinite(scenario.barrier_measure(iterate)));
   // the user-facing bounds of the problem (seen by the subproblem solver) are unchanged
   EXPECT_TRUE(std::isinf(scenario.barrier_problem.get_variables_lower_bounds()[0]));
   EXPECT_EQ(scenario.problem.get_variables_lower_bounds()[0], 0.);
}

TEST(BoundRelaxation, CenteredSlackFromMultiplier) {
   const double mu = 1e-2, z = 1e12; // mu/z = 1e-14 < cap
   BarrierScenario scenario(mu);
   Iterate iterate = scenario.make_iterate(0., 0.5, 1.5, z);
   Evaluations evaluations(scenario.model);
   scenario.barrier_problem.postprocess_iterate(iterate, evaluations);
   EXPECT_DOUBLE_EQ(scenario.barrier_problem.get_internal_variables_lower_bounds()[0], -mu / z);
}

TEST(BoundRelaxation, UpperBoundAndSlackVariable) {
   BarrierScenario scenario(0.1);
   // x1 slightly past its upper bound, s0 exactly on its lower bound (constraint bound 1)
   Iterate iterate = scenario.make_iterate(0.5, 1. + 1e-14, 1.);
   Evaluations evaluations(scenario.model);
   scenario.barrier_problem.postprocess_iterate(iterate, evaluations);
   EXPECT_EQ(scenario.barrier_problem.get_number_bound_relaxations(), 2u);
   EXPECT_EQ(iterate.primals[1], 1. + 1e-14);
   EXPECT_EQ(iterate.primals[2], 1.);
   EXPECT_GT(scenario.barrier_problem.get_internal_variables_upper_bounds()[1], iterate.primals[1]);
   EXPECT_LT(scenario.barrier_problem.get_internal_variables_lower_bounds()[2], 1.);
   EXPECT_TRUE(std::isfinite(scenario.barrier_measure(iterate)));
}

TEST(BoundRelaxation, RelaxationIsCumulativeAndMonotone) {
   BarrierScenario scenario(0.1);
   Evaluations evaluations(scenario.model);
   double previous_bound = 0.;
   for (size_t iteration = 0; iteration < 5; ++iteration) {
      // the iterate keeps landing on the current (relaxed) bound
      Iterate iterate = scenario.make_iterate(previous_bound, 0.5, 1.5);
      scenario.barrier_problem.postprocess_iterate(iterate, evaluations);
      const double bound = scenario.barrier_problem.get_internal_variables_lower_bounds()[0];
      EXPECT_LT(bound, previous_bound);
      EXPECT_GT(iterate.primals[0] - bound, 0.);
      previous_bound = bound;
   }
   EXPECT_EQ(scenario.barrier_problem.get_number_bound_relaxations(), 5u);
   // idempotent once the slack is safe
   Iterate iterate = scenario.make_iterate(previous_bound + slack_move, 0.5, 1.5);
   scenario.barrier_problem.postprocess_iterate(iterate, evaluations);
   EXPECT_EQ(scenario.barrier_problem.get_internal_variables_lower_bounds()[0], previous_bound);
}

TEST(BoundRelaxation, RelaxedBoundRemainsStrictlyBelowIterateDespiteRounding) {
   BarrierScenario scenario(1e-9);
   // s_min = eps*mu is below half an ulp of x = 1: x - s_min rounds to x without the nextafter safeguard
   Iterate iterate = scenario.make_iterate(0.5, 0.5, 1.);
   iterate.primals[2] = 1.; // on bound 1
   iterate.multipliers.lower_bounds[2] = 1e300;
   Evaluations evaluations(scenario.model);
   scenario.barrier_problem.postprocess_iterate(iterate, evaluations);
   EXPECT_LT(scenario.barrier_problem.get_internal_variables_lower_bounds()[2], iterate.primals[2]);
   EXPECT_TRUE(std::isfinite(scenario.barrier_measure(iterate)));
}

TEST(BoundRelaxation, PostprocessingDoesNotTriggerReevaluation) {
   BarrierScenario scenario(0.1);
   Iterate iterate = scenario.make_iterate(0., 0.5, 1.5);
   Evaluations evaluations(scenario.model);
   evaluations.evaluate_objective(scenario.model, iterate.primals);
   evaluations.evaluate_constraints(scenario.model, iterate.primals);
   const size_t objective_evaluations = scenario.model.number_objective_evaluations;
   const size_t constraint_evaluations = scenario.model.number_constraint_evaluations;
   scenario.barrier_problem.postprocess_iterate(iterate, evaluations);
   // the primals are unchanged: the cached evaluations remain valid
   evaluations.evaluate_objective(scenario.model, iterate.primals);
   evaluations.evaluate_constraints(scenario.model, iterate.primals);
   EXPECT_EQ(scenario.model.number_objective_evaluations, objective_evaluations);
   EXPECT_EQ(scenario.model.number_constraint_evaluations, constraint_evaluations);
}

TEST(BoundRelaxation, PhasesRelaxIndependently) {
   // two barrier problems over the same model, as in the two phases of FeasibilityRestoration
   BarrierScenario phase_a(0.1);
   const PrimalDualInteriorPointProblem phase_b(phase_a.problem, phase_a.parameters, phase_a.barrier_parameter, 100.);
   Evaluations evaluations(phase_a.model);

   // phase A relaxes its bound repeatedly; phase B is unaffected
   double x0 = 0.;
   for (size_t iteration = 0; iteration < 3; ++iteration) {
      Iterate iterate = phase_a.make_iterate(x0, 0.5, 1.5);
      phase_a.barrier_problem.postprocess_iterate(iterate, evaluations);
      x0 = phase_a.barrier_problem.get_internal_variables_lower_bounds()[0];
   }
   EXPECT_LT(x0, 0.);
   EXPECT_EQ(phase_b.get_internal_variables_lower_bounds()[0], 0.);

   // phase switch: the incoming iterate lies outside B's bounds; B relaxes just enough to contain it
   Iterate incoming = phase_a.make_iterate(x0, 0.5, 1.5);
   EXPECT_EQ(phase_b.relax_bounds_around(incoming), 1u);
   EXPECT_LT(phase_b.get_internal_variables_lower_bounds()[0], x0);
   EXPECT_GT(phase_b.get_internal_variables_lower_bounds()[0], x0 - 2. * slack_move);
   phase_b.set_auxiliary_measure(incoming);
   EXPECT_TRUE(std::isfinite(incoming.progress.auxiliary));
}

TEST(BoundRelaxation, SolveWithActiveBounds) {
   const BoundActiveModel model;
   Options options;
   DefaultOptions::load(options);
   Presets::set(options, "ipopt");
   options.set_double("primal_tolerance", 1e-12);
   options.set_double("dual_tolerance", 1e-12);
   Logger::level = SILENT;

   Uno uno;
   const Result result = uno.solve(model, options);
   ASSERT_EQ(result.optimization_status, OptimizationStatus::SUCCESS);
   EXPECT_NEAR(result.primal_solution[0], 0., 1e-8);
   EXPECT_NEAR(result.primal_solution[1], 1., 1e-8);
   // the solution may lie on a relaxed bound, but no further than the relaxation allows
   EXPECT_GE(result.primal_solution[0], -1e-8);
   EXPECT_LE(result.primal_solution[1], 1. + 1e-8);
}
