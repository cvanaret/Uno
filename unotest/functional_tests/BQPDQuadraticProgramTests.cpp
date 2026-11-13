// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "../DenseTestModel.hpp"
#include "ingredients/constraint_relaxation_strategies/relaxed_problems/l1RelaxedProblem.hpp"
#include "ingredients/hessian_models/ExactHessian.hpp"
#include "ingredients/inertia_correction_strategies/PrimalInertiaCorrection.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem_solvers/BQPD/BQPDQuadraticProgram.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

using namespace uno;

namespace {
   // min x0^2 + x1^2  s.t.  x0 + x1 <= 1.5, -10 <= x <= 10
   // (relaxed with one elastic variable p in the l1 feasibility problem)
   DenseTestModel convex_model() {
      return DenseTestModel("convex", 2, 1, ProblemType::NONLINEAR,
         [](const Vector<double>& x) { return x[0]*x[0] + x[1]*x[1]; },
         [](const Vector<double>& x, Vector<double>& g) { g[0] = 2.*x[0]; g[1] = 2.*x[1]; },
         [](const Vector<double>& x, Vector<double>& c) { c[0] = x[0] + x[1]; },
         [](const Vector<double>&, double* J) { J[0] = 1.; J[1] = 1.; },
         [](const Vector<double>&, double rho, const Vector<double>&, double* H) { H[0] = 2.*rho; H[1] = 0.; H[2] = 2.*rho; },
         {-10., -10.}, {10., 10.}, {-Inf}, {1.5}, {0.3, 0.4});
   }

   // vᵀ H v for a symmetric matrix stored as one COO triangle
   double explicit_quadratic_form(const BQPDQuadraticProgram& quadratic_program, const Vector<double>& vector) {
      double result = 0.;
      for (size_t nonzero_index = 0; nonzero_index < quadratic_program.hessian_values.size(); ++nonzero_index) {
         const size_t row = static_cast<size_t>(quadratic_program.hessian_row_indices[nonzero_index]);
         const size_t column = static_cast<size_t>(quadratic_program.hessian_column_indices[nonzero_index]);
         const double factor = (row == column) ? 1. : 2.;
         result += factor * quadratic_program.hessian_values[nonzero_index] * vector[row] * vector[column];
      }
      return result;
   }
} // namespace

// Bug: when the QP is solved with the explicit (regularized) Hessian, BQPDQuadraticProgram::compute_hessian_quadratic_form
// switched to the model's Hessian operator whenever the model provides one (e.g. AMPL models provide both). The
// quadratic form used for the predicted reduction was then not that of the QP that was solved: the operator adds the
// primal regularization to *all* the variables, including the elastic variables of the l1 feasibility problem that the
// explicit Hessian does not regularize, and every predicted reduction costs an extra Hessian-vector product.
TEST(BQPDQuadraticProgram, QuadraticFormUsesTheHessianOfTheQP) {
   if (SymmetricIndefiniteLinearSolverFactory::available_solvers().empty()) {
      GTEST_SKIP() << "the primal inertia correction requires a linear solver";
   }
   Logger::level = SILENT;
   DenseTestModel model = convex_model();
   model.enable_hessian_operator();
   Options options;
   DefaultOptions::load(options);

   // feasibility problem: zero objective multiplier and linear constraints => zero Hessian => regularization required
   const l1RelaxedProblem problem(model, 0., 1., true, false, 100.);
   ExactHessian hessian_model(model);
   PrimalInertiaCorrection inertia_correction_strategy(options);
   const Subproblem subproblem(problem, hessian_model, inertia_correction_strategy);
   ASSERT_EQ(subproblem.number_variables, 3); // x0, x1, p

   BQPDQuadraticProgram quadratic_program;
   quadratic_program.initialize_memory(subproblem);
   Iterate iterate(subproblem.number_variables, subproblem.number_constraints);
   iterate.primals[0] = 0.3;
   iterate.primals[1] = 0.4;
   iterate.primals[2] = 0.;
   Evaluations evaluations(model);
   WarmstartInformation warmstart_information{};
   Statistics statistics;
   inertia_correction_strategy.initialize_statistics(statistics);
   quadratic_program.fill(statistics, subproblem, iterate, /* trust_region_radius = */ Inf, evaluations, warmstart_information);
   ASSERT_GT(inertia_correction_strategy.get_primal_regularization_factor(), 0.);

   const Vector<double> vector{1., 1., 1.};
   const size_t number_hessian_vector_products = model.number_model_hessian_vector_products();
   const double quadratic_form = quadratic_program.compute_hessian_quadratic_form(subproblem, iterate, vector);
   EXPECT_DOUBLE_EQ(quadratic_form, explicit_quadratic_form(quadratic_program, vector));
   EXPECT_EQ(model.number_model_hessian_vector_products(), number_hessian_vector_products);
}