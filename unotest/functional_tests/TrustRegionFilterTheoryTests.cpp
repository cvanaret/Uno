// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Regression tests of the trust-region filter-SQP method against the global convergence theory of Fletcher, Leyffer &
// Toint, "On the global convergence of a filter-SQP algorithm", SIAM J. Optim. 13(1), 2002 (FLT)

#include <cmath>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
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
// non-KKT point, infinitely often, although every KKT point of the QP for radii in FLT's window [mu h, c sqrt(tau)]
// (tau: smallest infeasibility of the filter) is f-type and acceptable.
// Rule: if h <= theta (small-infeasibility threshold), an h-type iterate is accepted only if ||d||_inf <= b =
// kappa_h tau^nu with nu in (1/2, 1). Since h < tau, mu h << tau^nu << sqrt(tau) as tau -> 0: the radii [b, 2b] lie
// inside the window and the contradiction of Theorem 7 carries over to KKT points. After such a rejection, the radius
// decreases geometrically with b as a floor (max(b, min(radius, ||d||)/2)): jumping straight to b would skip the window
// used in the other case of Theorem 7 (finitely many h-type iterations).
// Here (default kappa_h = 100, nu = 0.75, theta = 1e-4): min -x - 1e-6 z  s.t.  x = 0, from (1e-5, 0). The QP step
// restores feasibility (the objective increases by 1e-5) and moves z to the trust-region bound (the objective decreases
// by 1e-6 radius): Delta q < 0 for every radius below 10, h-type step with ||d||_inf = radius. The initial upper bound
// of the filter (funnel width) is set just above h, so that tau = max(1e-5, factor 1e-5) is known and b ~ 0.02. From
// radius 1.5, the h-type steps of length 1.5, 0.75, ... > b are rejected, and the h-type step z = b is accepted. From
// radius 0.01 < b, the h-type step is accepted at once
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

   Options preset_options(const std::string& preset) {
      Options options;
      DefaultOptions::load(options);
      Presets::set(options, preset);
      options.set_string("QP_solver", "BQPD");
      options.set_string("LP_solver", "BQPD");
      Logger::level = SILENT;
      return options;
   }

   // solve one iteration from the radius; return (z, number of QPs)
   std::pair<double, size_t> solve_one_iteration(const DenseTestModel& model, Options& options, double radius) {
      options.set_double("TR_radius", radius);
      options.set_integer("max_iterations", 1);
      Uno uno;
      const Result result = uno.solve(model, options);
      EXPECT_EQ(result.number_iterations, 1) << "radius " << radius;
      EXPECT_NEAR(result.primal_solution[0], 0., 1e-12) << "radius " << radius;
      return {result.primal_solution[1], result.number_subproblems_solved};
   }

   void check_long_h_type_step_is_shortened(const std::string& preset) {
      const double h = 1e-5;
      const DenseTestModel model = h_type_step_model(h, 1e-6);
      for (const double initial_radius: {1.5, 0.01}) {
         Options options = preset_options(preset);
         // initial upper bound of the filter and width of the funnel: tau = max(ubd, factor h)
         options.set_double("filter_ubd", h);
         options.set_double("funnel_ubd", h);
         const double factor = options.get_double((preset == "funnelsqp") ? "funnel_fact" : "filter_fact");
         const double tau = std::max(h, factor*h);
         const double bound = 100.*std::pow(tau, 0.75);
         // expected radii: halving with the bound as a floor until the step (of length radius) is short enough
         size_t expected_number_QPs = 1;
         for (double radius = initial_radius; bound < radius; ++expected_number_QPs) {
            radius = std::max(bound, radius/2.);
         }
         const auto [z, number_QPs] = solve_one_iteration(model, options, initial_radius);
         EXPECT_NEAR(z, std::min(initial_radius, bound), 1e-12) << preset << ", radius " << initial_radius;
         EXPECT_EQ(number_QPs, expected_number_QPs) << preset << ", radius " << initial_radius;
      }
   }

   // above the small-infeasibility threshold, the bound does not apply. Here: min -x - 0.01 z  s.t.  x = 0, from
   // (0.1, 0) with radius 1.5: h = 0.1 > 1e-4, the h-type step z = 1.5 is accepted at once
   void check_h_type_step_above_threshold_is_not_shortened(const std::string& preset) {
      const DenseTestModel model = h_type_step_model(0.1, 0.01);
      Options options = preset_options(preset);
      const auto [z, number_QPs] = solve_one_iteration(model, options, 1.5);
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
// (h > small_infeasibility_threshold) as h-type. That is not an FLT h-type step (Delta q > 0), and it must not be
// shortened. The step-length bound only applies below the small-infeasibility threshold theta, which for the Waechter
// filter is its own small-infeasibility threshold: the threshold alone excludes such steps (the restriction of the bound
// to steps violating the switching condition is then implied).
// Here: min x - 0.01 z  s.t.  x = 0, from (0.1, 0) with radius 1.5 (trust region with the Waechter filter): Delta q =
// 0.1 + 0.015 > 0 satisfies the switching condition, and h = 0.1 > theta. The step z = 1.5 is accepted
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

#ifdef HAS_HIGHS
// Bug (MINLPTests nlp_cvx_expr_501_011, filtersqp with the convex QP solver HiGHS): with the bound b = kappa_h h^nu, the
// Newton steps of the end game were rejected as too long h-type steps. Near a KKT point on a curved constraint, moving by
// s along the constraint changes the infeasibility by O(s^2): the steps have length ~ sqrt(h) >> h^nu. The radius was
// cut, and HiGHS failed on the singular PSD Hessian. FLT's window is [mu h, c sqrt(tau)], so b = kappa_h tau^nu is
// enough for the theory; in a quadratically convergent end game, tau is the infeasibility of the previous h-type iterate
// ~ sqrt(h), and steps of length ~ sqrt(h) are accepted.
// Here: min -sum(x)  s.t.  sqrt(sum(x^2)) <= 1, from x = 1/n (n = 1, ..., 20). Solution x = 1/sqrt(n)
TEST(TrustRegionFilterTheory, NdSphereWithConvexQPSolver) {
   for (size_t n = 1; n <= 20; ++n) {
      const DenseTestModel model("nd_sphere", n, 1, ProblemType::NONLINEAR,
         [=](const Vector<double>& x) { double sum = 0.; for (size_t i = 0; i < n; ++i) { sum -= x[i]; } return sum; },
         [=](const Vector<double>&, Vector<double>& g) { for (size_t i = 0; i < n; ++i) { g[i] = -1.; } },
         [=](const Vector<double>& x, Vector<double>& c) {
            double sum = 0.; for (size_t i = 0; i < n; ++i) { sum += x[i]*x[i]; } c[0] = std::sqrt(sum); },
         [=](const Vector<double>& x, double* J) {
            double sum = 0.; for (size_t i = 0; i < n; ++i) { sum += x[i]*x[i]; }
            const double norm = std::sqrt(sum); for (size_t i = 0; i < n; ++i) { J[i] = x[i]/norm; } },
         // Hessian of the Lagrangian rho f - y c (lower triangle, row by row): -y (I/r - x x^T/r^3)
         [=](const Vector<double>& x, double, const Vector<double>& y, double* H) {
            double sum = 0.; for (size_t i = 0; i < n; ++i) { sum += x[i]*x[i]; }
            const double norm = std::sqrt(sum);
            size_t index = 0;
            for (size_t i = 0; i < n; ++i) {
               for (size_t j = 0; j <= i; ++j) {
                  H[index++] = -y[0]*(((i == j) ? 1./norm : 0.) - x[i]*x[j]/(norm*norm*norm));
               }
            }
         },
         std::vector<double>(n, -Inf), std::vector<double>(n, Inf), {-Inf}, {1.}, std::vector<double>(n, 1./static_cast<double>(n)));
      Options options = filtersqp_options();
      options.set_string("QP_solver", "HiGHS");
      options.set_integer("max_iterations", 10000);
      Uno uno;
      const Result result = uno.solve(model, options);
      EXPECT_EQ(result.optimization_status, OptimizationStatus::SUCCESS) << "n = " << n;
      // |f(x) - f*| <= ||grad f|| dist(x, x*) ~ sqrt(n) primal_tolerance
      const double tolerance = std::sqrt(static_cast<double>(n))*options.get_double("primal_tolerance");
      EXPECT_NEAR(result.solution_objective, -std::sqrt(static_cast<double>(n)), tolerance) << "n = " << n;
   }
}
#endif

// Bug: at an infeasible iterate of the optimality phase, a radius below TR_min_radius threw "Small radius" (algorithmic
// error). In FLT's algorithm, reducing the radius at an infeasible point always ends with an incompatible QP (Lemma 6),
// that is, an h-type iteration via the restoration phase, never with a failure. In floating-point arithmetic, the QP
// may remain compatible for tiny radii (e.g. when the feasibility tolerance of the QP solver accepts a tiny box). The
// method should switch to the feasibility problem instead.
// Here: min x  s.t.  x^2 = 1, from x = 0.1 with radius 10 and a large TR_min_radius (3) that triggers the small-radius
// logic deterministically: the linearized constraint requires d = 4.95, the trial iterate 5.05 increases both the
// objective and the infeasibility and is rejected, and the radius becomes 4.95/2 < 3
TEST(TrustRegionFilterTheory, SmallRadiusAtInfeasibleIterateSwitchesToRestoration) {
   const DenseTestModel model("small_radius_infeasible", 1, 1, ProblemType::NONLINEAR,
      [](const Vector<double>& x) { return x[0]; },
      [](const Vector<double>&, Vector<double>& g) { g[0] = 1.; },
      [](const Vector<double>& x, Vector<double>& c) { c[0] = x[0]*x[0]; },
      [](const Vector<double>& x, double* J) { J[0] = 2.*x[0]; },
      [](const Vector<double>&, double, const Vector<double>& y, double* H) { H[0] = -2.*y[0]; },
      {-Inf}, {Inf}, {1.}, {1.}, {0.1});
   Options options = filtersqp_options();
   options.set_double("TR_radius", 10.);
   options.set_double("TR_min_radius", 3.);

   Uno uno;
   const Result result = uno.solve(model, options);
   ASSERT_EQ(result.optimization_status, OptimizationStatus::SUCCESS);
   EXPECT_LE(result.solution_primal_feasibility, 1e-6);
   EXPECT_NEAR(std::abs(result.primal_solution[0]), 1., 1e-6);
}
