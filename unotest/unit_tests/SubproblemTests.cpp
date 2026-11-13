// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "../DenseTestModel.hpp"
#include "ingredients/constraint_relaxation_strategies/relaxed_problems/l1RelaxedProblem.hpp"
#include "ingredients/hessian_models/ExactHessian.hpp"
#include "ingredients/inertia_correction_strategies/PrimalInertiaCorrection.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "optimization/Iterate.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

using namespace uno;

// Bug: Subproblem::compute_hessian_vector_product added the primal regularization to *all* the variables, whereas the
// explicit regularized Hessian (compute_regularized_hessian_sparsity/regularize_lagrangian_hessian) only regularizes the
// primal regularization variables. In the l1 feasibility problem, the elastic variables were regularized by the operator
// but not by the matrix.
TEST(Subproblem, HessianVectorProductRegularizesTheSameVariablesAsTheMatrix) {
   if (SymmetricIndefiniteLinearSolverFactory::available_solvers().empty()) {
      GTEST_SKIP() << "the primal inertia correction requires a linear solver";
   }
   Logger::level = SILENT;
   // min x0^2 + x1^2  s.t.  x0 + x1 <= 1.5
   DenseTestModel model("convex", 2, 1, ProblemType::NONLINEAR,
      [](const Vector<double>& x) { return x[0]*x[0] + x[1]*x[1]; },
      [](const Vector<double>& x, Vector<double>& g) { g[0] = 2.*x[0]; g[1] = 2.*x[1]; },
      [](const Vector<double>& x, Vector<double>& c) { c[0] = x[0] + x[1]; },
      [](const Vector<double>&, double* J) { J[0] = 1.; J[1] = 1.; },
      [](const Vector<double>&, double rho, const Vector<double>&, double* H) { H[0] = 2.*rho; H[1] = 0.; H[2] = 2.*rho; },
      {-10., -10.}, {10., 10.}, {-Inf}, {1.5}, {0.3, 0.4});
   model.enable_hessian_operator();
   Options options;
   DefaultOptions::load(options);

   // feasibility problem (x0, x1, p): zero objective multiplier and linear constraint => zero Hessian, regularized
   const l1RelaxedProblem problem(model, 0., 1., true, false, 100.);
   ExactHessian hessian_model(model);
   PrimalInertiaCorrection inertia_correction_strategy(options);
   const Subproblem subproblem(problem, hessian_model, inertia_correction_strategy);
   ASSERT_EQ(subproblem.number_variables, 3);
   ASSERT_EQ(subproblem.get_primal_regularization_variables().size(), 2); // the elastic variable is not regularized

   Iterate iterate(subproblem.number_variables, subproblem.number_constraints);
   iterate.primals[0] = 0.3;
   iterate.primals[1] = 0.4;
   Statistics statistics;
   inertia_correction_strategy.initialize_statistics(statistics);
   Vector<double> hessian_values(subproblem.number_regularized_hessian_nonzeros());
   subproblem.evaluate_lagrangian_hessian(statistics, iterate, hessian_values.view());
   subproblem.regularize_lagrangian_hessian(statistics, hessian_values.view());
   const double regularization_factor = inertia_correction_strategy.get_primal_regularization_factor();
   ASSERT_GT(regularization_factor, 0.);

   const Vector<double> vector{1., 1., 1.};
   Vector<double> result(subproblem.number_variables);
   subproblem.compute_hessian_vector_product(iterate, iterate.primals.view(), vector.view(), result.view());
   EXPECT_DOUBLE_EQ(result[0], regularization_factor);
   EXPECT_DOUBLE_EQ(result[1], regularization_factor);
   EXPECT_DOUBLE_EQ(result[2], 0.); // elastic variable: no curvature, no regularization
}