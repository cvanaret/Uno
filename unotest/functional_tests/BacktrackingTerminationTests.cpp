// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// At a feasible iterate (h = 0), the minimum step length of Wächter & Biegler (eq. (21)) is 0, and the backtracking line
// search must eventually accept a step (Remark 1): the Armijo condition holds for a small enough step length. In floating
// point, once alpha*d is below the resolution of x, the trial iterates coincide with the current iterate, but the
// acceptance test still depends on alpha through the predicted reduction (roundoff guards, roundoff-feasible steps).
// Backtracking must therefore go on, but without evaluating the functions again at the same point. Without the roundoff
// guards (protect_actual_reduction_against_roundoff, roundoff_feasibility_factor) and without an absolute Armijo
// tolerance (which accepts the null step once the predicted reduction is below it), the line search used to evaluate the
// functions at every halving until the step length underflowed (about 1075 trial iterates per iteration)

#include <cmath>
#include <gtest/gtest.h>
#include "Uno.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "optimization/Result.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"
#include "tools/Logger.hpp"
#include "../DenseTestModel.hpp"

using namespace uno;

namespace {
   // min x0^2 + x1^2 + 1 s.t. x0 - x1 = 0, from (1, 1). The linear constraint remains satisfied exactly (h = 0), and
   // the user gradient has the wrong sign: the Newton direction is an ascent direction and the Armijo condition never
   // holds. Returns the number of objective evaluations
   size_t number_objective_evaluations_with_wrong_gradient(bool roundoff_guards) {
      size_t number_evaluations = 0;
      const DenseTestModel model("wrong_gradient", 2, 1, ProblemType::NONLINEAR,
         [&](const Vector<double>& x) { ++number_evaluations; return x[0] * x[0] + x[1] * x[1] + 1.; },
         [](const Vector<double>& x, Vector<double>& gradient) { gradient[0] = -2. * x[0]; gradient[1] = -2. * x[1]; },
         [](const Vector<double>& x, Vector<double>& constraints) { constraints[0] = x[0] - x[1]; },
         [](const Vector<double>&, double* jacobian) { jacobian[0] = 1.; jacobian[1] = -1.; },
         [](const Vector<double>&, double objective_multiplier, const Vector<double>&, double* hessian) {
            hessian[0] = 2. * objective_multiplier; hessian[1] = 0.; hessian[2] = 2. * objective_multiplier;
         },
         {-INFINITY, -INFINITY}, {INFINITY, INFINITY}, {0.}, {0.}, {1., 1.}, {0});
      Options options;
      DefaultOptions::load(options);
      Presets::set(options, "ipopt");
      options.set_bool("protect_actual_reduction_against_roundoff", roundoff_guards);
      options.set_double("roundoff_feasibility_factor", roundoff_guards ? 1e-2 : 0.);
      options.set_double("armijo_tolerance", 0.);
      options.set_integer("max_iterations", 50);
      Logger::level = SILENT;
      Uno uno;
      const Result result = uno.solve(model, options);
      EXPECT_NE(result.solution_status, SolutionStatus::FEASIBLE_KKT_POINT);
      return number_evaluations;
   }
} // namespace

TEST(BacktrackingTermination, CoincidingTrialIteratesAreNotReevaluated) {
   if (SymmetricIndefiniteLinearSolverFactory::available_solvers().empty()) {
      GTEST_SKIP() << "no linear solver available";
   }
   // about 55 halvings bring alpha*d below the resolution of x = 1: the functions are evaluated about 55 times per iteration
   EXPECT_LT(number_objective_evaluations_with_wrong_gradient(false), 300u);
}

namespace {
   // F(t) = 1e-3 + G s + H s^2/2 + M s^4 with s = t - 1: at t = 1, the Newton step is d = -G/H = -1e-14, a descent
   // direction just above the resolution of x = 1. The quartic term rejects the full step: backtracking brings alpha*d
   // below the resolution of x, where the trial iterates coincide with the current iterate
   constexpr double G = 1e7;
   constexpr double H = 1e21;
   constexpr double M = 1e57;
   double F(double t) { const double s = t - 1.; return 1e-3 + G * s + 0.5 * H * s * s + M * s * s * s * s; }
   double dF(double t) { const double s = t - 1.; return G + H * s + 4. * M * s * s * s; }
   double d2F(double t) { const double s = t - 1.; return H + 12. * M * s * s; }
} // namespace

// a step is accepted once alpha is small enough (roundoff guard of the Armijo condition): the line search must not stop
// when the trial iterate coincides with the current iterate. Same outcome as without the reuse of the evaluations
TEST(BacktrackingTermination, DescentDirectionBelowResolutionIsEventuallyAccepted) {
   if (SymmetricIndefiniteLinearSolverFactory::available_solvers().empty()) {
      GTEST_SKIP() << "no linear solver available";
   }
   const DenseTestModel model("quartic", 2, 1, ProblemType::NONLINEAR,
      [](const Vector<double>& x) { return F(x[0]) + F(x[1]); },
      [](const Vector<double>& x, Vector<double>& gradient) { gradient[0] = dF(x[0]); gradient[1] = dF(x[1]); },
      [](const Vector<double>& x, Vector<double>& constraints) { constraints[0] = x[0] - x[1]; },
      [](const Vector<double>&, double* jacobian) { jacobian[0] = 1.; jacobian[1] = -1.; },
      [](const Vector<double>& x, double objective_multiplier, const Vector<double>&, double* hessian) {
         hessian[0] = objective_multiplier * d2F(x[0]); hessian[1] = 0.; hessian[2] = objective_multiplier * d2F(x[1]);
      },
      {-INFINITY, -INFINITY}, {INFINITY, INFINITY}, {0.}, {0.}, {1., 1.}, {0});
   Options options;
   DefaultOptions::load(options);
   Presets::set(options, "ipopt");
   options.set_integer("max_iterations", 6);
   Logger::level = SILENT;
   Uno uno;
   const Result result = uno.solve(model, options);
   EXPECT_EQ(result.optimization_status, OptimizationStatus::SUCCESS);
   EXPECT_EQ(result.solution_status, SolutionStatus::FEASIBLE_SMALL_STEP);
}

// with the roundoff guards (IPOPT preset), the number of evaluations is also bounded
TEST(BacktrackingTermination, RoundoffGuardsBoundTheNumberOfEvaluations) {
   if (SymmetricIndefiniteLinearSolverFactory::available_solvers().empty()) {
      GTEST_SKIP() << "no linear solver available";
   }
   EXPECT_LT(number_objective_evaluations_with_wrong_gradient(true), 300u);
}
