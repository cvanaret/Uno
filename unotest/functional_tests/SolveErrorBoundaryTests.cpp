// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Regression tests: errors raised before the optimization loop (initial point, function scaling) must be reported in
// the Result, not escape Uno::solve()

#include <cmath>
#include <stdexcept>
#include <gtest/gtest.h>
#include "../DenseTestModel.hpp"
#include "Uno.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "optimization/Result.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"
#include "tools/Logger.hpp"

using namespace uno;

namespace {
   enum class Failure {INITIAL_POINT_RUNTIME_ERROR, INITIAL_POINT_EVALUATION_ERROR, GRADIENT_RUNTIME_ERROR};

   // min (x0 - 1)^2 + (x1 - 2)^2  s.t.  x0 + x1 <= 1
   class FailingModel: public DenseTestModel {
   public:
      explicit FailingModel(Failure failure):
         DenseTestModel("failing_model", 2, 1, ProblemType::NONLINEAR,
            [](const Vector<double>& x) { return std::pow(x[0] - 1., 2) + std::pow(x[1] - 2., 2); },
            [](const Vector<double>& x, Vector<double>& g) { g[0] = 2.*(x[0] - 1.); g[1] = 2.*(x[1] - 2.); },
            [](const Vector<double>& x, Vector<double>& c) { c[0] = x[0] + x[1]; },
            [](const Vector<double>&, double* J) { J[0] = 1.; J[1] = 1.; },
            [](const Vector<double>&, double rho, const Vector<double>&, double* H) { H[0] = 2.*rho; H[1] = 0.; H[2] = 2.*rho; },
            {-Inf, -Inf}, {Inf, Inf}, {-Inf}, {1.}, {0., 0.}),
         failure(failure) { }

      void initial_primal_point(Vector<double>& x) const override {
         if (this->failure == Failure::INITIAL_POINT_RUNTIME_ERROR) {
            throw std::runtime_error("no initial point");
         }
         if (this->failure == Failure::INITIAL_POINT_EVALUATION_ERROR) {
            throw FunctionEvaluationError();
         }
         DenseTestModel::initial_primal_point(x);
      }

      // a user callback that throws a std::exception that is not an EvaluationError
      void evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const override {
         if (this->failure == Failure::GRADIENT_RUNTIME_ERROR) {
            throw std::runtime_error("user callback failure");
         }
         DenseTestModel::evaluate_objective_gradient(x, gradient);
      }

   private:
      const Failure failure;
   };

   struct Outcome {
      bool threw{true};
      OptimizationStatus status{OptimizationStatus::SUCCESS};
      double objective{0.};
   };

   Outcome solve(const Model& model, const Options& options) {
      Outcome outcome;
      try {
         Uno uno;
         const Result result = uno.solve(model, options);
         outcome = {false, result.optimization_status, result.solution_objective};
      }
      catch (...) {
      }
      return outcome;
   }

   Options default_options() {
      Options options;
      DefaultOptions::load(options);
      Presets::set(options, "ipopt");
      Logger::level = SILENT;
      return options;
   }
} // namespace

TEST(SolveErrorBoundary, InitialPointRuntimeError) {
   const FailingModel model(Failure::INITIAL_POINT_RUNTIME_ERROR);
   const Outcome outcome = solve(model, default_options());
   ASSERT_FALSE(outcome.threw);
   EXPECT_EQ(outcome.status, OptimizationStatus::ALGORITHMIC_ERROR);
   EXPECT_TRUE(std::isnan(outcome.objective));
}

TEST(SolveErrorBoundary, InitialPointEvaluationError) {
   const FailingModel model(Failure::INITIAL_POINT_EVALUATION_ERROR);
   const Outcome outcome = solve(model, default_options());
   ASSERT_FALSE(outcome.threw);
   EXPECT_EQ(outcome.status, OptimizationStatus::EVALUATION_ERROR);
}

// Bug: ScaledModel caught only EvaluationError: another exception thrown by a user callback during the scaling escaped
// Uno::solve()
TEST(SolveErrorBoundary, NonEvaluationErrorDuringScaling) {
   const FailingModel model(Failure::GRADIENT_RUNTIME_ERROR);
   Options options = default_options();
   options.set_bool("use_function_scaling", true);
   const Outcome outcome = solve(model, options);
   ASSERT_FALSE(outcome.threw);
   EXPECT_EQ(outcome.status, OptimizationStatus::ALGORITHMIC_ERROR);
}
