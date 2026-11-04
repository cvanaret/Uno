// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Regression tests: a failure of the final evaluation of the user functions (after the postprocessing of the solution)
// must not be reported as a success, and the Result must not carry stale or partially written function values.

#include <cmath>
#include <gtest/gtest.h>
#include "../DenseTestModel.hpp"
#include "Uno.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Result.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"
#include "tools/Logger.hpp"

using namespace uno;

namespace {
   // min (x0 - 1)^2 + (x1 - 2)^2  s.t.  x0 + x1 <= 1
   // The postprocessing invalidates the cached evaluations (as the snapping of fixed variables does): the functions
   // are evaluated again at the final iterate, where the selected ones fail.
   class FailingFinalEvaluationModel: public DenseTestModel {
   public:
      FailingFinalEvaluationModel(bool objective_fails, bool constraints_fail):
         DenseTestModel("failing_final_evaluation", 2, 1, ProblemType::NONLINEAR,
            [](const Vector<double>& x) { return std::pow(x[0] - 1., 2) + std::pow(x[1] - 2., 2); },
            [](const Vector<double>& x, Vector<double>& g) { g[0] = 2.*(x[0] - 1.); g[1] = 2.*(x[1] - 2.); },
            [](const Vector<double>& x, Vector<double>& c) { c[0] = x[0] + x[1]; },
            [](const Vector<double>&, double* J) { J[0] = 1.; J[1] = 1.; },
            [](const Vector<double>&, double rho, const Vector<double>&, double* H) { H[0] = 2.*rho; H[1] = 0.; H[2] = 2.*rho; },
            {-Inf, -Inf}, {Inf, Inf}, {-Inf}, {1.}, {0., 0.}),
         objective_fails(objective_fails), constraints_fail(constraints_fail) { }

      [[nodiscard]] double evaluate_objective(const Vector<double>& x) const override {
         if (this->postprocessed && this->objective_fails) {
            throw FunctionEvaluationError();
         }
         return DenseTestModel::evaluate_objective(x);
      }

      void evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const override {
         if (this->postprocessed && this->constraints_fail) {
            constraints[0] = 42.; // partially written
            throw FunctionEvaluationError();
         }
         DenseTestModel::evaluate_constraints(x, constraints);
      }

      void postprocess_solution(Iterate& /*iterate*/, Evaluations& evaluations) const override {
         evaluations.reset();
         this->postprocessed = true;
      }

   private:
      const bool objective_fails;
      const bool constraints_fail;
      mutable bool postprocessed{false};
   };

   Options ipopt_options() {
      Options options;
      DefaultOptions::load(options);
      Presets::set(options, "ipopt");
      Logger::level = SILENT;
      return options;
   }
} // namespace

TEST(FinalEvaluation, ValidFinalEvaluationIsSuccess) {
   if (SymmetricIndefiniteLinearSolverFactory::available_solvers().empty()) {
      GTEST_SKIP() << "the interior-point method requires a linear solver";
   }
   const FailingFinalEvaluationModel model(false, false);
   Uno uno;
   const Result result = uno.solve(model, ipopt_options());
   ASSERT_EQ(result.optimization_status, OptimizationStatus::SUCCESS);
   EXPECT_NEAR(result.solution_objective, 2., 1e-6); // optimum (0, 1)
   EXPECT_NEAR(result.constraint_values[0], 1., 1e-6);
}

TEST(FinalEvaluation, ObjectiveFailureIsEvaluationError) {
   if (SymmetricIndefiniteLinearSolverFactory::available_solvers().empty()) {
      GTEST_SKIP() << "the interior-point method requires a linear solver";
   }
   const FailingFinalEvaluationModel model(true, false);
   Uno uno;
   const Result result = uno.solve(model, ipopt_options());
   EXPECT_EQ(result.optimization_status, OptimizationStatus::EVALUATION_ERROR);
   EXPECT_TRUE(std::isnan(result.solution_objective));
   EXPECT_NEAR(result.constraint_values[0], 1., 1e-6);
}

TEST(FinalEvaluation, ConstraintFailureIsEvaluationError) {
   if (SymmetricIndefiniteLinearSolverFactory::available_solvers().empty()) {
      GTEST_SKIP() << "the interior-point method requires a linear solver";
   }
   const FailingFinalEvaluationModel model(false, true);
   Uno uno;
   const Result result = uno.solve(model, ipopt_options());
   EXPECT_EQ(result.optimization_status, OptimizationStatus::EVALUATION_ERROR);
   EXPECT_NEAR(result.solution_objective, 2., 1e-6);
   // no partially written value
   EXPECT_TRUE(std::isnan(result.constraint_values[0]));
}
