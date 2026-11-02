// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// The feasibility restoration rescales the constraints of the feasibility problem where it starts: a ScaledModel with
// identity scaling wraps the (scaled) model of the original problem, and its constraints are rescaled on demand. The
// unscaled quantities must be expressed in the units of the user model through both layers.

#include <gtest/gtest.h>
#include "../DenseTestModel.hpp"
#include "linear_algebra/Norm.hpp"
#include "model/ScaledModel.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"

using namespace uno;

namespace {
   // 1e4 x0 + x1 = 2e4 and x0 - x1 <= 1
   DenseTestModel badly_scaled_model() {
      return DenseTestModel("badly_scaled", 2, 2, ProblemType::NONLINEAR,
         [](const Vector<double>& x) { return x[0]*x[0] + x[1]*x[1]; },
         [](const Vector<double>& x, Vector<double>& g) { g[0] = 2.*x[0]; g[1] = 2.*x[1]; },
         [](const Vector<double>& x, Vector<double>& c) { c[0] = 1e4*x[0] + x[1]; c[1] = x[0] - x[1]; },
         [](const Vector<double>&, double* J) { J[0] = 1e4; J[1] = 1.; J[2] = 1.; J[3] = -1.; },
         [](const Vector<double>&, double rho, const Vector<double>&, double* H) { H[0] = 2.*rho; H[1] = 0.; H[2] = 2.*rho; },
         {-Inf, -Inf}, {Inf, Inf}, {2e4, -Inf}, {2e4, 1.}, {0., 0.});
   }
} // namespace

TEST(ScaledModel, RescaledConstraintsAndNestedUnscaling) {
   Options options;
   DefaultOptions::load(options);
   const DenseTestModel model = badly_scaled_model();
   ScaledModel scaled_model(model, options);
   // identity scaling until the constraints are rescaled
   EXPECT_EQ(scaled_model.get_constraint_scaling()[0], 1.);
   EXPECT_EQ(scaled_model.get_constraints_upper_bounds()[0], 2e4);

   const Vector<double> jacobian_values{1e4, 1., 1., -1.};
   scaled_model.rescale_constraints(jacobian_values);
   // threshold 100 on the row norms of the Jacobian
   EXPECT_DOUBLE_EQ(scaled_model.get_constraint_scaling()[0], 1e-2);
   EXPECT_DOUBLE_EQ(scaled_model.get_constraint_scaling()[1], 1.);
   EXPECT_DOUBLE_EQ(scaled_model.get_constraints_lower_bounds()[0], 200.);
   EXPECT_DOUBLE_EQ(scaled_model.get_constraints_upper_bounds()[0], 200.);

   // nested scaling: the outer model scales the rows of the inner (scaled) model by 1/2
   ScaledModel nested_model(scaled_model, options);
   const Vector<double> scaled_jacobian_values{100., 1e-2, 200., -200.};
   nested_model.rescale_constraints(scaled_jacobian_values);
   EXPECT_DOUBLE_EQ(nested_model.get_constraint_scaling()[0], 1.);
   EXPECT_DOUBLE_EQ(nested_model.get_constraint_scaling()[1], 0.5);

   // at x = (1, 0): c(x) = (1e4, 1), violation of the user constraints (1e4, 0)
   Vector<double> x{1., 0.};
   Vector<double> constraints(2);
   nested_model.evaluate_constraints(x, constraints);
   EXPECT_DOUBLE_EQ(constraints[0], 100.);
   EXPECT_DOUBLE_EQ(constraints[1], 0.5);
   EXPECT_DOUBLE_EQ(nested_model.unscaled_primal_infeasibility(0., constraints, Norm::INF), 1e4);
   // the bounds are scaled through both layers
   EXPECT_DOUBLE_EQ(nested_model.get_constraints_upper_bounds()[1], 0.5);
}
