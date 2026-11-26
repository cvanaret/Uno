// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Tests of the user callbacks of the C API: contract of the return codes, and what the callbacks receive

#include <gtest/gtest.h>
#include <string>
#include <vector>
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

      uno_model* model{nullptr};
      uno_solver* solver{nullptr};
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
   uno_model* create_failing_hs015_model(EvaluationState& state) {
      double variables_lower_bounds[] = {-Inf, -Inf};
      double variables_upper_bounds[] = {0.5, Inf};
      double constraints_lower_bounds[] = {1., 0.};
      double constraints_upper_bounds[] = {Inf, Inf};
      uno_int jacobian_row_indices[] = {0, 1, 0, 1};
      uno_int jacobian_column_indices[] = {0, 0, 1, 1};
      uno_int hessian_row_indices[] = {0, 1, 1};
      uno_int hessian_column_indices[] = {0, 0, 1};
      double x0[] = {-2., 1.};
      uno_model* model = uno_create_model(UNO_PROBLEM_NONLINEAR, number_variables, variables_lower_bounds,
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
         uno_model* model = create_failing_hs015_model(state);
         uno_int status{-1};
         EXPECT_EQ(uno_optimize(this->solver, model), UNO_OPTIMIZE_SUCCESS);
         EXPECT_EQ(uno_get_optimization_status(this->solver, &status), UNO_RESULT_SUCCESS);
         EXPECT_GE(state.number_calls, 2) << "the function was not evaluated twice";
         uno_destroy_model(model);
         return status;
      }

      uno_solver* solver{nullptr};
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

// the solver callbacks receive the iterate in the space of the user model: "number_variables" variables and
// "number_constraints" constraints (no slacks or elastics), and multipliers with the sign convention of the user
namespace {
   struct RecordedIterate {
      uno_int number_variables{-1};
      uno_int number_constraints{-1};
      std::vector<double> primals;
      std::vector<double> lower_bound_multipliers;
      std::vector<double> upper_bound_multipliers;
      std::vector<double> constraint_multipliers;
   };

   struct IterateRecorder {
      std::vector<RecordedIterate> notified;
      std::vector<RecordedIterate> termination_checks;
   };

   RecordedIterate record(uno_int number_variables, uno_int number_constraints, const double* primals,
         const double* lower_bound_multipliers, const double* upper_bound_multipliers, const double* constraint_multipliers) {
      const auto n = static_cast<size_t>(number_variables);
      const auto m = static_cast<size_t>(number_constraints);
      return {number_variables, number_constraints, {primals, primals + n},
         {lower_bound_multipliers, lower_bound_multipliers + n}, {upper_bound_multipliers, upper_bound_multipliers + n},
         {constraint_multipliers, constraint_multipliers + m}};
   }

   void recording_notify(uno_int number_variables, uno_int number_constraints, const double* primals,
         const double* lower_bound_multipliers, const double* upper_bound_multipliers, const double* constraint_multipliers,
         double /*objective_multiplier*/, double /*primal_feasibility_residual*/, double /*stationarity_residual*/,
         double /*complementarity_residual*/, void* user_data) {
      static_cast<IterateRecorder*>(user_data)->notified.push_back(record(number_variables, number_constraints, primals,
         lower_bound_multipliers, upper_bound_multipliers, constraint_multipliers));
   }

   uno_int recording_termination(uno_int number_variables, uno_int number_constraints, const double* primals,
         const double* lower_bound_multipliers, const double* upper_bound_multipliers, const double* constraint_multipliers,
         double /*objective_multiplier*/, double /*primal_feasibility_residual*/, double /*stationarity_residual*/,
         double /*complementarity_residual*/, void* user_data) {
      static_cast<IterateRecorder*>(user_data)->termination_checks.push_back(record(number_variables, number_constraints,
         primals, lower_bound_multipliers, upper_bound_multipliers, constraint_multipliers));
      return 0;
   }

   // min x0 s.t. x0^2 + x1^2 <= -1: infeasible, Uno ends in the feasibility restoration phase
   uno_int infeasible_objective(uno_int /*n*/, const double* x, double* objective_value, void* /*user_data*/) {
      *objective_value = x[0];
      return 0;
   }

   uno_int infeasible_gradient(uno_int /*n*/, const double* /*x*/, double* gradient, void* /*user_data*/) {
      gradient[0] = 1.;
      gradient[1] = 0.;
      return 0;
   }

   uno_int infeasible_constraints(uno_int /*n*/, uno_int /*m*/, const double* x, double* values, void* /*user_data*/) {
      values[0] = x[0]*x[0] + x[1]*x[1];
      return 0;
   }

   uno_int infeasible_jacobian(uno_int /*n*/, uno_int /*nnz*/, const double* x, double* values, void* /*user_data*/) {
      values[0] = 2.*x[0];
      values[1] = 2.*x[1];
      return 0;
   }

   uno_int infeasible_hessian(uno_int /*n*/, uno_int /*m*/, uno_int /*nnz*/, const double* /*x*/,
         double /*objective_multiplier*/, const double* multipliers, double* values, void* /*user_data*/) {
      // Lagrangian rho f(x) - y c(x) (negative sign convention)
      values[0] = -2.*multipliers[0];
      values[1] = -2.*multipliers[0];
      return 0;
   }

   uno_model* create_infeasible_model() {
      double variables_lower_bounds[] = {-10., -10.};
      double variables_upper_bounds[] = {10., 10.};
      double constraints_lower_bounds[] = {-Inf};
      double constraints_upper_bounds[] = {-1.};
      uno_int jacobian_row_indices[] = {0, 0};
      uno_int jacobian_column_indices[] = {0, 1};
      uno_int hessian_indices[] = {0, 1};
      double x0[] = {1., 2.};
      uno_model* model = uno_create_model(UNO_PROBLEM_NONLINEAR, 2, variables_lower_bounds, variables_upper_bounds,
         UNO_ZERO_BASED_INDEXING);
      uno_set_objective(model, UNO_MINIMIZE, infeasible_objective, infeasible_gradient);
      uno_set_constraints(model, 1, infeasible_constraints, constraints_lower_bounds, constraints_upper_bounds, 2,
         jacobian_row_indices, jacobian_column_indices, infeasible_jacobian);
      uno_set_lagrangian_hessian(model, 2, UNO_LOWER_TRIANGLE, hessian_indices, hessian_indices, infeasible_hessian);
      uno_set_initial_primal_iterate(model, x0);
      return model;
   }

   void expect_user_dimensions(const std::vector<RecordedIterate>& iterates, uno_int n, uno_int m) {
      ASSERT_FALSE(iterates.empty());
      for (const RecordedIterate& iterate: iterates) {
         EXPECT_EQ(iterate.number_variables, n);
         EXPECT_EQ(iterate.number_constraints, m);
      }
   }

   class CAPICallbackIterateTests: public ::testing::TestWithParam<uno_int> {
   protected:
      void SetUp() override {
         if (available_preset() == nullptr) {
            GTEST_SKIP() << "no subproblem solver available";
         }
         this->solver = uno_create_solver();
         ASSERT_NE(this->solver, nullptr);
         ASSERT_EQ(uno_set_solver_string_option(this->solver, "logger", "SILENT"), UNO_OPTIONS_SUCCESS);
         ASSERT_EQ(uno_set_solver_preset(this->solver, available_preset()), UNO_OPTIONS_SUCCESS);
         uno_set_solver_callbacks(this->solver, recording_notify, recording_termination, &this->recorder);
      }

      void TearDown() override {
         uno_destroy_solver(this->solver);
      }

      uno_solver* solver{nullptr};
      IterateRecorder recorder;
   };
} // namespace

// hs015 has 2 variables and 2 inequality constraints (slacked by the interior-point method)
TEST_P(CAPICallbackIterateTests, CallbacksReceiveUserDimensions) {
   uno_model* model = create_hs015_model();
   uno_set_lagrangian_sign_convention(model, GetParam());
   ASSERT_EQ(uno_optimize(this->solver, model), UNO_OPTIMIZE_SUCCESS);
   uno_destroy_model(model);
   expect_user_dimensions(this->recorder.notified, number_variables, number_constraints);
   expect_user_dimensions(this->recorder.termination_checks, number_variables, number_constraints);
}

// at convergence, the last accepted iterate is the solution: the notify callback sees the same values as the result getters
TEST_P(CAPICallbackIterateTests, LastNotifiedIterateIsTheSolution) {
   uno_model* model = create_hs015_model();
   uno_set_lagrangian_sign_convention(model, GetParam());
   ASSERT_EQ(uno_optimize(this->solver, model), UNO_OPTIMIZE_SUCCESS);
   uno_destroy_model(model);
   uno_int status{-1};
   ASSERT_EQ(uno_get_optimization_status(this->solver, &status), UNO_RESULT_SUCCESS);
   ASSERT_EQ(status, UNO_SUCCESS);
   ASSERT_FALSE(this->recorder.notified.empty());
   const RecordedIterate& last = this->recorder.notified.back();
   ASSERT_EQ(last.number_variables, number_variables);
   ASSERT_EQ(last.number_constraints, number_constraints);

   std::vector<double> primals(2), lower_bound_duals(2), upper_bound_duals(2), constraint_duals(2);
   ASSERT_EQ(uno_get_primal_solution(this->solver, primals.data()), UNO_RESULT_SUCCESS);
   ASSERT_EQ(uno_get_lower_bound_dual_solution(this->solver, lower_bound_duals.data()), UNO_RESULT_SUCCESS);
   ASSERT_EQ(uno_get_upper_bound_dual_solution(this->solver, upper_bound_duals.data()), UNO_RESULT_SUCCESS);
   ASSERT_EQ(uno_get_constraint_dual_solution(this->solver, constraint_duals.data()), UNO_RESULT_SUCCESS);
   for (size_t index: {0, 1}) {
      EXPECT_DOUBLE_EQ(last.primals[index], primals[index]);
      EXPECT_DOUBLE_EQ(last.lower_bound_multipliers[index], lower_bound_duals[index]);
      EXPECT_DOUBLE_EQ(last.upper_bound_multipliers[index], upper_bound_duals[index]);
      EXPECT_DOUBLE_EQ(last.constraint_multipliers[index], constraint_duals[index]);
   }
}

// in the feasibility restoration phase, the elastic variables are not exposed
TEST_P(CAPICallbackIterateTests, CallbacksReceiveUserDimensionsDuringRestoration) {
   uno_model* model = create_infeasible_model();
   uno_set_lagrangian_sign_convention(model, GetParam());
   ASSERT_EQ(uno_optimize(this->solver, model), UNO_OPTIMIZE_SUCCESS);
   uno_destroy_model(model);
   expect_user_dimensions(this->recorder.notified, 2, 1);
   expect_user_dimensions(this->recorder.termination_checks, 2, 1);
}

namespace {
   std::string sign_convention_name(const ::testing::TestParamInfo<uno_int>& info) {
      return (info.param == UNO_MULTIPLIER_NEGATIVE) ? "negative" : "positive";
   }
} // namespace

INSTANTIATE_TEST_SUITE_P(SignConventions, CAPICallbackIterateTests, ::testing::Values(UNO_MULTIPLIER_NEGATIVE,
   UNO_MULTIPLIER_POSITIVE), sign_convention_name);
