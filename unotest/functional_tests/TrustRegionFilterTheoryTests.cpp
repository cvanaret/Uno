// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Regression tests of the trust-region filter-SQP method against the global convergence theory of Fletcher, Leyffer &
// Toint, "On the global convergence of a filter-SQP algorithm", SIAM J. Optim. 13(1), 2002 (FLT)

#include <cmath>
#include <string>
#include <tuple>
#include <utility>
#include <gtest/gtest.h>
#include "../DenseTestModel.hpp"
#include "Uno.hpp"
#include "optimization/Result.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"
#include "tools/Logger.hpp"

using namespace uno;

namespace {
   Options filtersqp_options() {
      Options options;
      DefaultOptions::load(options);
      Presets::set(options, "filtersqp");
      options.set_string("QP_solver", "BQPD");
      options.set_string("LP_solver", "BQPD");
      Logger::level = SILENT;
      return options;
   }
} // namespace

// Bug: an h-type iterate was accepted whatever the step length. FLT's proof (Theorem 7, h-type case) relies on the inner
// loop sampling an f-type step before an h-type step: Delta q decreases monotonically as the radius decreases, which is
// a consequence of the *global* optimality of the QP solution. BQPD returns KKT points of nonconvex QPs: at a large
// radius, a KKT point with Delta q <= delta h^2 could be accepted as an h-type iterate in the neighborhood of a feasible
// non-KKT point, infinitely often, although every KKT point of the QP for radii in [mu h, kappa] is f-type and acceptable.
// Rule: if h <= theta (small-infeasibility threshold), an h-type iterate is accepted only if ||d||_inf <= b =
// kappa_h h^nu with nu in (1/2, 1). Since mu h << h^nu << sqrt(h) for small h, the radii [b, 2b] lie inside FLT's window
// and the contradiction of Theorem 7 carries over to KKT points; along the subsequence of Theorem 7, h -> 0 drops below
// any fixed theta > 0. After such a rejection, the radius decreases geometrically with b as a floor
// (max(b, min(radius, ||d||)/2)): jumping straight to b would skip the window [mu h, rho_bar] used in the other case of
// Theorem 7 (finitely many h-type iterations), which requires the accepted radii to be bounded away from 0.
// Here (default kappa_h = 100, nu = 0.75, theta = 1e-4): min -x - 1e-6 z  s.t.  x = 0, from (1e-5, 0). The QP step
// restores feasibility (the objective increases by 1e-5) and moves z to the trust-region bound (the objective decreases
// by 1e-6 radius): Delta q < 0 for every radius below 10, h-type step with ||d||_inf = radius, and b = 100 (1e-5)^0.75
// ~ 0.0178. From radius 1.5, the h-type steps of length 1.5, 0.75, ..., 0.0234375 are rejected, and the h-type step
// z = b is accepted (8 QPs). From radius 0.01 <= b, the h-type step is accepted at once
namespace {
   DenseTestModel h_type_step_model(double x0, double z_weight) {
      return DenseTestModel("h_type_step", 2, 1, ProblemType::NONLINEAR,
         [=](const Vector<double>& x) { return -x[0] - z_weight*x[1]; },
         [=](const Vector<double>&, Vector<double>& g) { g[0] = -1.; g[1] = -z_weight; },
         [](const Vector<double>& x, Vector<double>& c) { c[0] = x[0]; },
         [](const Vector<double>&, double* J) { J[0] = 1.; J[1] = 0.; },
         [](const Vector<double>&, double, const Vector<double>&, double* H) { H[0] = 0.; H[1] = 0.; H[2] = 0.; },
         {-Inf, -10.}, {Inf, 10.}, {0.}, {0.}, {x0, 0.}, {0});
   }

   // solve one iteration from the radius; return (z, number of QPs)
   std::pair<double, size_t> solve_one_iteration(const DenseTestModel& model, const std::string& preset, double radius) {
      Options options;
      DefaultOptions::load(options);
      Presets::set(options, preset);
      options.set_string("QP_solver", "BQPD");
      options.set_string("LP_solver", "BQPD");
      Logger::level = SILENT;
      options.set_double("TR_radius", radius);
      options.set_integer("max_iterations", 1);

      Uno uno;
      const Result result = uno.solve(model, options);
      EXPECT_EQ(result.number_iterations, 1) << preset << ", radius " << radius;
      EXPECT_NEAR(result.primal_solution[0], 0., 1e-12) << preset << ", radius " << radius;
      return {result.primal_solution[1], result.number_subproblems_solved};
   }

   void check_long_h_type_step_is_shortened(const std::string& preset) {
      const DenseTestModel model = h_type_step_model(1e-5, 1e-6);
      const double bound = 100.*std::pow(1e-5, 0.75);
      for (const auto& [radius, number_subproblems]: {std::pair{1.5, size_t(8)}, std::pair{0.01, size_t(1)}}) {
         const auto [z, number_QPs] = solve_one_iteration(model, preset, radius);
         EXPECT_NEAR(z, std::min(radius, bound), 1e-12) << preset << ", radius " << radius;
         // geometric decrease of the radius with the bound as a floor
         EXPECT_EQ(number_QPs, number_subproblems) << preset << ", radius " << radius;
      }
   }

   // above the small-infeasibility threshold, the bound does not apply. Here: min -x - 0.01 z  s.t.  x = 0, from
   // (0.1, 0) with radius 1.5: h = 0.1 > 1e-4, the h-type step z = 1.5 is accepted at once
   void check_h_type_step_above_threshold_is_not_shortened(const std::string& preset) {
      const DenseTestModel model = h_type_step_model(0.1, 0.01);
      const auto [z, number_QPs] = solve_one_iteration(model, preset, 1.5);
      EXPECT_NEAR(z, 1.5, 1e-12) << preset;
      EXPECT_EQ(number_QPs, 1) << preset;
   }
} // namespace

TEST(TrustRegionFilterTheory, LongHTypeStepIsShortened) {
   check_long_h_type_step_is_shortened("filtersqp");
}

TEST(TrustRegionFilterTheory, LongHTypeStepIsShortenedFunnel) {
   check_long_h_type_step_is_shortened("funnelsqp");
}

TEST(TrustRegionFilterTheory, HTypeStepAboveThresholdIsNotShortened) {
   check_h_type_step_above_threshold_is_not_shortened("filtersqp");
   check_h_type_step_above_threshold_is_not_shortened("funnelsqp");
}

// The Waechter filter treats a trial iterate that satisfies the switching condition at a non-small infeasibility
// (h > small_infeasibility_threshold) as h-type. That is not an FLT h-type step (Delta q > 0), and FLT's argument only
// concerns infeasibilities h -> 0, where the infeasibility is small: the step-length bound only applies if the
// switching condition is violated.
// Here: min x - 0.01 z  s.t.  x = 0, from (0.1, 0) with radius 1.5 (trust region with the Waechter filter): Delta q =
// 0.1 + 0.015 > 0 satisfies the switching condition and h = 0.1 is not small. The step z = 1.5 is accepted
TEST(TrustRegionFilterTheory, WaechterFilterLongSwitchingStepIsNotShortened) {
   const DenseTestModel model("long_switching_step", 2, 1, ProblemType::NONLINEAR,
      [](const Vector<double>& x) { return x[0] - 0.01*x[1]; },
      [](const Vector<double>&, Vector<double>& g) { g[0] = 1.; g[1] = -0.01; },
      [](const Vector<double>& x, Vector<double>& c) { c[0] = x[0]; },
      [](const Vector<double>&, double* J) { J[0] = 1.; J[1] = 0.; },
      [](const Vector<double>&, double, const Vector<double>&, double* H) { H[0] = 0.; H[1] = 0.; H[2] = 0.; },
      {-Inf, -10.}, {Inf, 10.}, {0.}, {0.}, {0.1, 0.}, {0});
   Options options = filtersqp_options();
   options.set_string("globalization_strategy", "waechter_filter_method");
   options.set_double("TR_radius", 1.5);
   options.set_integer("max_iterations", 1);

   Uno uno;
   const Result result = uno.solve(model, options);
   ASSERT_EQ(result.number_iterations, 1);
   EXPECT_NEAR(result.primal_solution[0], 0., 1e-12);
   EXPECT_NEAR(result.primal_solution[1], 1.5, 1e-12);
   EXPECT_EQ(result.number_subproblems_solved, 1);
}

// the exponent of the h-type step bound must lie in (1/2, 1), the factor must be positive, and the small-infeasibility
// threshold below which the bound applies must be positive (otherwise, the bound never applies and case (i) of
// Theorem 7 is open again)
TEST(TrustRegionFilterTheory, HTypeStepBoundOptionsAreValidated) {
   const DenseTestModel model = h_type_step_model(1e-5, 1e-6);
   for (const auto& [factor, exponent, threshold, valid]: {std::tuple{100., 0.75, 1e-4, true}, std::tuple{100., 0.5, 1e-4, false},
         std::tuple{100., 1., 1e-4, false}, std::tuple{0., 0.75, 1e-4, false}, std::tuple{100., 0.75, 0., false}}) {
      Options options = filtersqp_options();
      options.set_double("TR_h_type_step_factor", factor);
      options.set_double("TR_h_type_step_exponent", exponent);
      options.set_double("barrier_small_infeasibility_factor", threshold);
      options.set_integer("max_iterations", 1);
      Uno uno;
      const Result result = uno.solve(model, options);
      EXPECT_EQ(result.optimization_status != OptimizationStatus::ALGORITHMIC_ERROR, valid) << factor << ", " << exponent << ", " << threshold;
   }
}
