// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Regression tests of the trust-region filter-SQP method against the global convergence theory of Fletcher, Leyffer &
// Toint, "On the global convergence of a filter-SQP algorithm", SIAM J. Optim. 13(1), 2002 (FLT)

#include <cmath>
#include <string>
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
// Rule: an h-type iterate is accepted only if ||d||_inf <= kappa_h h^nu with nu in (1/2, 1). Since mu h << h^nu << sqrt(h)
// for small h, the radii [kappa_h h^nu, 2 kappa_h h^nu] lie inside FLT's window and the contradiction of Theorem 7
// carries over to KKT points. Otherwise, the radius is set to kappa_h h^nu and the QP is solved again.
// Here (default kappa_h = 1, nu = 0.75): min -x - 0.01 z  s.t.  x = 0, from (0.1, 0) with radius 1.5. The QP step restores feasibility (the objective
// increases by 0.1) and moves z to the trust-region bound (the objective decreases by 0.015): Delta q < 0, h-type step
// with ||d||_inf = 1.5 > 0.1^0.75. The step is rejected, the radius becomes 0.1^0.75 and the h-type step z = 0.1^0.75
// is accepted
TEST(TrustRegionFilterTheory, LongHTypeStepIsShortened) {
   const DenseTestModel model("long_h_type_step", 2, 1, ProblemType::NONLINEAR,
      [](const Vector<double>& x) { return -x[0] - 0.01*x[1]; },
      [](const Vector<double>&, Vector<double>& g) { g[0] = -1.; g[1] = -0.01; },
      [](const Vector<double>& x, Vector<double>& c) { c[0] = x[0]; },
      [](const Vector<double>&, double* J) { J[0] = 1.; J[1] = 0.; },
      [](const Vector<double>&, double, const Vector<double>&, double* H) { H[0] = 0.; H[1] = 0.; H[2] = 0.; },
      {-Inf, -10.}, {Inf, 10.}, {0.}, {0.}, {0.1, 0.}, {0});
   for (const double radius: {1.5, 0.1}) {
      Options options = filtersqp_options();
      options.set_double("TR_radius", radius);
      options.set_integer("max_iterations", 1);

      Uno uno;
      const Result result = uno.solve(model, options);
      ASSERT_EQ(result.number_iterations, 1) << radius;
      EXPECT_NEAR(result.primal_solution[0], 0., 1e-12) << radius;
      // a short h-type step (radius 0.1 <= 0.1^0.75) is accepted at once
      EXPECT_NEAR(result.primal_solution[1], std::min(radius, std::pow(0.1, 0.75)), 1e-12) << radius;
   }
}
