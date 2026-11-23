// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "../DenseTestModel.hpp"
#include "ingredients/constraint_relaxation_strategies/relaxed_problems/l1RelaxedProblem.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"

using namespace uno;

// Bug: the predicted reduction of the objective of the l1 feasibility problem (sum of the elastic variables) ignored the
// curvature term -1/2 d^T H d, whereas the restoration QP minimizes rho sum(p + n) + 1/2 d^T H d, with H the Hessian of
// the Lagrangian of the feasibility problem (sum of multipliers times constraint Hessians, generally indefinite). The
// switching condition and the Armijo test (FLT: Delta q = q(0) - q(d)) were then applied to the decrease of a different
// model: a nonconvex restoration QP step with a positive Delta q could be classified as h-type (linear part <= 0).
TEST(l1RelaxedProblemPredictedReduction, IncludesCurvature) {
   // x0^2 + x1 <= 1 (one elastic variable p: x0^2 + x1 - p <= 1)
   const DenseTestModel model("feasibility", 2, 1, ProblemType::NONLINEAR,
      [](const Vector<double>&) { return 0.; },
      [](const Vector<double>&, Vector<double>& g) { g[0] = 0.; g[1] = 0.; },
      [](const Vector<double>& x, Vector<double>& c) { c[0] = x[0]*x[0] + x[1]; },
      [](const Vector<double>& x, double* J) { J[0] = 2.*x[0]; J[1] = 1.; },
      [](const Vector<double>&, double, const Vector<double>& y, double* H) { H[0] = -2.*y[0]; H[1] = 0.; H[2] = 0.; },
      {-10., -10.}, {10., 10.}, {-Inf}, {1.}, {1., 1.});
   const double coefficient = 2.;
   const l1RelaxedProblem problem(model, 0., coefficient, true, false, 100.);
   ASSERT_EQ(problem.number_variables, 3);

   const Iterate iterate(problem.number_variables, problem.number_constraints);
   Evaluations evaluations(model);
   // the elastic variable decreases by 0.5
   const Vector<double> direction{0.1, -0.2, -0.5};
   const double hessian_quadratic_form = 3.;
   const auto predicted_objective_reduction = problem.build_predicted_objective_reduction(iterate, direction,
      evaluations, hessian_quadratic_form);
   const double linear_reduction = coefficient*0.5;
   for (const double step_length: {1., 0.5}) {
      EXPECT_DOUBLE_EQ(predicted_objective_reduction(step_length, 0., false), step_length*linear_reduction);
      EXPECT_DOUBLE_EQ(predicted_objective_reduction(step_length, 0., true),
         step_length*linear_reduction - step_length*step_length/2.*hessian_quadratic_form);
   }
}
