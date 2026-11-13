// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// A user callback stops Uno immediately by throwing UserTermination

#include <cmath>
#include <string>
#include <vector>
#include <gtest/gtest.h>
#include "../DenseTestModel.hpp"
#include "Uno.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "optimization/Multipliers.hpp"
#include "optimization/Result.hpp"
#include "optimization/UserTermination.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"
#include "tools/Logger.hpp"
#include "tools/UserCallbacks.hpp"

using namespace uno;

namespace {
   enum class Failure { USER_TERMINATION, EVALUATION_ERROR };

   // hs015: min 100 (x1 - x0^2)^2 + (1 - x0)^2  s.t.  x0 x1 >= 1, x0 + x1^2 >= 0, x0 <= 0.5. The objective throws
   // at its failing_call-th evaluation (1-based, 0 never fails): a UserTermination from then on, a FunctionEvaluationError
   // only once
   class InterruptibleModel: public DenseTestModel {
   public:
      InterruptibleModel(): DenseTestModel("hs015", 2, 2, ProblemType::NONLINEAR,
            [](const Vector<double>& x) { return 100.*(x[1] - x[0]*x[0])*(x[1] - x[0]*x[0]) + (1. - x[0])*(1. - x[0]); },
            [](const Vector<double>& x, Vector<double>& g) {
               g[0] = 400.*x[0]*x[0]*x[0] - 400.*x[0]*x[1] + 2.*x[0] - 2.;
               g[1] = 200.*(x[1] - x[0]*x[0]);
            },
            [](const Vector<double>& x, Vector<double>& c) { c[0] = x[0]*x[1]; c[1] = x[0] + x[1]*x[1]; },
            [](const Vector<double>& x, double* J) { J[0] = x[1]; J[1] = x[0]; J[2] = 1.; J[3] = 2.*x[1]; },
            [](const Vector<double>& x, double rho, const Vector<double>& y, double* H) {
               H[0] = rho*(1200.*x[0]*x[0] - 400.*x[1] + 2.);
               H[1] = -400.*rho*x[0] - y[0];
               H[2] = 200.*rho - 2.*y[1];
            },
            {-Inf, -Inf}, {0.5, Inf}, {1., 0.}, {Inf, Inf}, {-2., 1.}) {
      }

      [[nodiscard]] double evaluate_objective(const Vector<double>& x) const override {
         ++this->number_calls;
         if (0 < this->failing_call) {
            if (this->failure == Failure::USER_TERMINATION && this->failing_call <= this->number_calls) {
               throw UserTermination();
            }
            if (this->failure == Failure::EVALUATION_ERROR && this->failing_call == this->number_calls) {
               throw FunctionEvaluationError();
            }
         }
         return DenseTestModel::evaluate_objective(x);
      }

      Failure failure{Failure::USER_TERMINATION};
      size_t failing_call{0};
      mutable size_t number_calls{0};
   };

   // records the accepted iterates (and the number of objective evaluations of the model when they are accepted); the
   // termination callback is never true
   class RecordingCallbacks: public UserCallbacks {
   public:
      explicit RecordingCallbacks(const InterruptibleModel& model): model(model) { }

      void notify_acceptable_iterate(const Vector<double>& primals, const Multipliers& /*multipliers*/,
            double /*objective_multiplier*/, double /*primal_feasibility_residual*/, double /*stationarity_residual*/,
            double /*complementarity_residual*/) override {
         this->accepted_primals.emplace_back(primals[0], primals[1]);
         this->number_objective_calls.push_back(this->model.number_calls);
         if (this->throw_at_notification == this->accepted_primals.size()) {
            throw UserTermination();
         }
      }
      bool termination(const Vector<double>& /*primals*/, const Multipliers& /*multipliers*/, double /*objective_multiplier*/,
            double /*primal_feasibility_residual*/, double /*stationarity_residual*/, double /*complementarity_residual*/) override {
         ++this->number_termination_calls;
         return false;
      }

      const InterruptibleModel& model;
      std::vector<std::pair<double, double>> accepted_primals{};
      std::vector<size_t> number_objective_calls{};
      size_t throw_at_notification{0}; // 0: never
      size_t number_termination_calls{0};
   };

   Options preset_options(const std::string& preset) {
      Options options;
      DefaultOptions::load(options);
      Presets::set(options, preset);
      Logger::level = SILENT;
      return options;
   }

   std::vector<std::string> available_presets() {
      std::vector<std::string> presets;
      if (!SymmetricIndefiniteLinearSolverFactory::available_solvers().empty()) {
         presets.emplace_back("ipopt");
      }
#ifdef HAS_BQPD
      presets.emplace_back("filtersqp");
#endif
      return presets;
   }
} // namespace

// Feature: a model evaluation that throws UserTermination stops Uno at once (no backtracking, no further evaluation)
// with USER_TERMINATION, at the last accepted iterate
TEST(UserTermination, ThrownByModelStopsImmediately) {
   const std::vector<std::string> presets = available_presets();
   if (presets.empty()) {
      GTEST_SKIP() << "no preset available";
   }
   for (const std::string& preset: presets) {
      // reference run: number of objective evaluations when each iterate is accepted
      InterruptibleModel reference_model;
      RecordingCallbacks reference_callbacks(reference_model);
      Uno reference_uno;
      const Result reference = reference_uno.solve(reference_model, preset_options(preset), reference_callbacks);
      ASSERT_EQ(reference.optimization_status, OptimizationStatus::SUCCESS) << preset;
      ASSERT_GE(reference_callbacks.accepted_primals.size(), 3) << preset;

      // the first objective evaluation after the second accepted iterate requests the termination
      InterruptibleModel model;
      model.failing_call = reference_callbacks.number_objective_calls[1] + 1;
      RecordingCallbacks callbacks(model);
      Uno uno;
      const Result result = uno.solve(model, preset_options(preset), callbacks);
      EXPECT_EQ(result.optimization_status, OptimizationStatus::USER_TERMINATION) << preset;
      EXPECT_EQ(result.number_iterations, 3) << preset; // the third iteration was interrupted
      // the result is the last accepted iterate
      ASSERT_EQ(callbacks.accepted_primals.size(), 2) << preset;
      EXPECT_EQ(result.primal_solution[0], callbacks.accepted_primals.back().first) << preset;
      EXPECT_EQ(result.primal_solution[1], callbacks.accepted_primals.back().second) << preset;
      EXPECT_EQ(callbacks.accepted_primals.back(), reference_callbacks.accepted_primals[1]) << preset;
      // the objective of the last accepted iterate is reported (it was computed before the termination request)
      EXPECT_TRUE(std::isfinite(result.solution_objective)) << preset;
      EXPECT_EQ(model.number_calls, model.failing_call) << preset;
   }
}

// contrast: an EvaluationError at the same evaluation is handled by the globalization mechanism, which backtracks
TEST(UserTermination, EvaluationErrorDoesNotStopUno) {
   const std::vector<std::string> presets = available_presets();
   if (presets.empty()) {
      GTEST_SKIP() << "no preset available";
   }
   for (const std::string& preset: presets) {
      InterruptibleModel reference_model;
      RecordingCallbacks reference_callbacks(reference_model);
      Uno reference_uno;
      static_cast<void>(reference_uno.solve(reference_model, preset_options(preset), reference_callbacks));
      ASSERT_GE(reference_callbacks.accepted_primals.size(), 3) << preset;

      // the objective fails once, at the same evaluation as above
      InterruptibleModel model;
      model.failure = Failure::EVALUATION_ERROR;
      model.failing_call = reference_callbacks.number_objective_calls[1] + 1;
      Uno uno;
      const Result result = uno.solve(model, preset_options(preset));
      EXPECT_EQ(result.optimization_status, OptimizationStatus::SUCCESS) << preset;
      EXPECT_GT(result.number_iterations, 3) << preset;
      EXPECT_GT(model.number_calls, model.failing_call + 1) << preset;
   }
}

// a termination request at the initial point is reported as USER_TERMINATION, not as an evaluation error
TEST(UserTermination, ThrownAtInitialPoint) {
   for (const std::string& preset: available_presets()) {
      InterruptibleModel model;
      model.failing_call = 1;
      Uno uno;
      const Result result = uno.solve(model, preset_options(preset));
      EXPECT_EQ(result.optimization_status, OptimizationStatus::USER_TERMINATION) << preset;
      EXPECT_EQ(result.number_iterations, 0) << preset;
      // the objective is not called again (e.g. to evaluate the final iterate)
      EXPECT_EQ(model.number_calls, 1) << preset;
   }
}

// the solver callbacks can also stop Uno immediately
TEST(UserTermination, ThrownByNotifyCallback) {
   for (const std::string& preset: available_presets()) {
      InterruptibleModel model;
      RecordingCallbacks callbacks(model);
      callbacks.throw_at_notification = 2;
      Uno uno;
      const Result result = uno.solve(model, preset_options(preset), callbacks);
      EXPECT_EQ(result.optimization_status, OptimizationStatus::USER_TERMINATION) << preset;
      EXPECT_EQ(callbacks.accepted_primals.size(), 2) << preset;
      // the iterate that was being notified was not accepted yet: the result is the previously accepted iterate
      EXPECT_EQ(result.primal_solution[0], callbacks.accepted_primals.front().first) << preset;
      EXPECT_EQ(result.primal_solution[1], callbacks.accepted_primals.front().second) << preset;
   }
}
