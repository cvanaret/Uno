// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Tests of the user callbacks of the C API: contract of the return codes, and what the callbacks receive

#include <gtest/gtest.h>
#include <string>
#include "CAPITestProblem.hpp"

using namespace uno;
using namespace uno::capi_test;

namespace {
   // state shared with the termination callback through user_data
   struct TerminationState {
      uno_int return_value{0};
      int number_calls{0};
   };

   uno_int termination_callback(uno_int /*number_variables*/, uno_int /*number_constraints*/, const double* /*primals*/,
         const double* /*lower_bound_multipliers*/, const double* /*upper_bound_multipliers*/,
         const double* /*constraint_multipliers*/, double /*objective_multiplier*/, double /*primal_feasibility_residual*/,
         double /*stationarity_residual*/, double /*complementarity_residual*/, void* user_data) {
      auto* state = static_cast<TerminationState*>(user_data);
      ++state->number_calls;
      return state->return_value;
   }

   class CAPICallbackTests: public ::testing::Test {
   protected:
      void SetUp() override {
         if (available_preset() == nullptr) {
            GTEST_SKIP() << "no subproblem solver available";
         }
         this->model = create_hs015_model();
         this->solver = uno_create_solver();
         ASSERT_NE(this->model, nullptr);
         ASSERT_NE(this->solver, nullptr);
         ASSERT_EQ(uno_set_solver_string_option(this->solver, "logger", "SILENT"), UNO_OPTIONS_SUCCESS);
         ASSERT_EQ(uno_set_solver_preset(this->solver, available_preset()), UNO_OPTIONS_SUCCESS);
      }

      void TearDown() override {
         uno_destroy_solver(this->solver);
         uno_destroy_model(this->model);
      }

      [[nodiscard]] uno_int optimization_status() const {
         uno_int status{-1};
         EXPECT_EQ(uno_get_optimization_status(this->solver, &status), UNO_RESULT_SUCCESS);
         return status;
      }

      void* model{nullptr};
      void* solver{nullptr};
   };
} // namespace

// a termination callback that returns 0 lets the optimization run to convergence
TEST_F(CAPICallbackTests, TerminationCallbackReturningZeroContinues) {
   TerminationState state{0, 0};
   uno_set_solver_callbacks(this->solver, nullptr, termination_callback, &state);
   ASSERT_EQ(uno_optimize(this->solver, this->model), UNO_OPTIMIZE_SUCCESS);
   EXPECT_EQ(this->optimization_status(), UNO_SUCCESS);
   EXPECT_GT(state.number_calls, 1);
}

// a termination callback that returns a nonzero value stops the optimization at the first call
TEST_F(CAPICallbackTests, TerminationCallbackReturningNonzeroStops) {
   TerminationState state{1, 0};
   uno_set_solver_callbacks(this->solver, nullptr, termination_callback, &state);
   ASSERT_EQ(uno_optimize(this->solver, this->model), UNO_OPTIMIZE_SUCCESS);
   EXPECT_EQ(this->optimization_status(), UNO_USER_TERMINATION);
   EXPECT_EQ(state.number_calls, 1);
}

// return codes of the evaluation callbacks: 0 = success, positive = evaluation error, negative = user termination
namespace {
   enum class Function {OBJECTIVE, CONSTRAINTS, GRADIENT, JACOBIAN, HESSIAN};

   // state shared with the evaluation callbacks through the model user_data
   struct EvaluationState {
      Function failing_function{Function::OBJECTIVE};
      int failing_call{0}; // the failing_call-th evaluation of failing_function returns return_code
      uno_int return_code{0};
      int number_calls{0};
   };

   uno_int next_return_code(void* user_data, Function function) {
      auto* state = static_cast<EvaluationState*>(user_data);
      if (state->failing_function == function) {
         ++state->number_calls;
         if (state->number_calls == state->failing_call) {
            return state->return_code;
         }
      }
      return 0;
   }

   uno_int failing_objective(uno_int n, const double* x, double* objective_value, void* user_data) {
      objective_function(n, x, objective_value, user_data);
      return next_return_code(user_data, Function::OBJECTIVE);
   }

   uno_int failing_constraints(uno_int n, uno_int m, const double* x, double* values, void* user_data) {
      constraint_functions(n, m, x, values, user_data);
      return next_return_code(user_data, Function::CONSTRAINTS);
   }

   uno_int failing_gradient(uno_int n, const double* x, double* gradient, void* user_data) {
      objective_gradient(n, x, gradient, user_data);
      return next_return_code(user_data, Function::GRADIENT);
   }

   uno_int failing_jacobian(uno_int n, uno_int nnz, const double* x, double* values, void* user_data) {
      jacobian(n, nnz, x, values, user_data);
      return next_return_code(user_data, Function::JACOBIAN);
   }

   uno_int failing_hessian(uno_int n, uno_int m, uno_int nnz, const double* x, double objective_multiplier,
         const double* multipliers, double* values, void* user_data) {
      lagrangian_hessian(n, m, nnz, x, objective_multiplier, multipliers, values, user_data);
      return next_return_code(user_data, Function::HESSIAN);
   }

   // hs015 whose callbacks may fail, according to the EvaluationState passed as user_data
   void* create_failing_hs015_model(EvaluationState& state) {
      double variables_lower_bounds[] = {-Inf, -Inf};
      double variables_upper_bounds[] = {0.5, Inf};
      double constraints_lower_bounds[] = {1., 0.};
      double constraints_upper_bounds[] = {Inf, Inf};
      uno_int jacobian_row_indices[] = {0, 1, 0, 1};
      uno_int jacobian_column_indices[] = {0, 0, 1, 1};
      uno_int hessian_row_indices[] = {0, 1, 1};
      uno_int hessian_column_indices[] = {0, 0, 1};
      double x0[] = {-2., 1.};
      void* model = uno_create_model(UNO_PROBLEM_NONLINEAR, number_variables, variables_lower_bounds,
         variables_upper_bounds, UNO_ZERO_BASED_INDEXING);
      uno_set_objective(model, UNO_MINIMIZE, failing_objective, failing_gradient);
      uno_set_constraints(model, number_constraints, failing_constraints, constraints_lower_bounds,
         constraints_upper_bounds, 4, jacobian_row_indices, jacobian_column_indices, failing_jacobian);
      uno_set_lagrangian_hessian(model, 3, UNO_LOWER_TRIANGLE, hessian_row_indices, hessian_column_indices,
         failing_hessian);
      uno_set_initial_primal_iterate(model, x0);
      uno_set_user_data(model, &state);
      return model;
   }

   class CAPIEvaluationReturnCodeTests: public ::testing::TestWithParam<Function> {
   protected:
      void SetUp() override {
         if (available_preset() == nullptr) {
            GTEST_SKIP() << "no subproblem solver available";
         }
         this->solver = uno_create_solver();
         ASSERT_NE(this->solver, nullptr);
         ASSERT_EQ(uno_set_solver_string_option(this->solver, "logger", "SILENT"), UNO_OPTIONS_SUCCESS);
         ASSERT_EQ(uno_set_solver_preset(this->solver, available_preset()), UNO_OPTIONS_SUCCESS);
      }

      void TearDown() override {
         uno_destroy_solver(this->solver);
      }

      // solves hs015 where the second evaluation of the function GetParam() returns "return_code"
      [[nodiscard]] uno_int solve(uno_int return_code) {
         EvaluationState state{GetParam(), 2, return_code, 0};
         void* model = create_failing_hs015_model(state);
         uno_int status{-1};
         EXPECT_EQ(uno_optimize(this->solver, model), UNO_OPTIMIZE_SUCCESS);
         EXPECT_EQ(uno_get_optimization_status(this->solver, &status), UNO_RESULT_SUCCESS);
         EXPECT_GE(state.number_calls, 2) << "the function was not evaluated twice";
         uno_destroy_model(model);
         return status;
      }

      void* solver{nullptr};
   };
} // namespace

// a negative return code is a user termination request
TEST_P(CAPIEvaluationReturnCodeTests, NegativeReturnCodeRequestsTermination) {
   EXPECT_EQ(this->solve(-1), UNO_USER_TERMINATION);
}

// a zero return code is a success: the optimization converges
TEST_P(CAPIEvaluationReturnCodeTests, ZeroReturnCodeIsSuccess) {
   EXPECT_EQ(this->solve(0), UNO_SUCCESS);
}

// a positive return code is an evaluation error, never a user termination
TEST_P(CAPIEvaluationReturnCodeTests, PositiveReturnCodeIsNotATermination) {
   EXPECT_NE(this->solve(1), UNO_USER_TERMINATION);
}

namespace {
   std::string function_name(const ::testing::TestParamInfo<Function>& info) {
      constexpr const char* names[] = {"objective", "constraints", "gradient", "jacobian", "hessian"};
      return names[static_cast<int>(info.param)];
   }
} // namespace

INSTANTIATE_TEST_SUITE_P(AllEvaluations, CAPIEvaluationReturnCodeTests, ::testing::Values(Function::OBJECTIVE,
   Function::CONSTRAINTS, Function::GRADIENT, Function::JACOBIAN, Function::HESSIAN), function_name);
