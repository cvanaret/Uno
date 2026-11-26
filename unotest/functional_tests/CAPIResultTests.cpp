// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Tests of the functions of the C API that get the result of an optimization and of the model functions: they report
// errors with status codes (or false/nullptr) and never throw

#include <gtest/gtest.h>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>
#include "CAPITestProblem.hpp"
#include "tools/Logger.hpp"

using namespace uno;
using namespace uno::capi_test;

// the declarations are noexcept for C++ callers
static_assert(noexcept(uno_get_version(nullptr, nullptr, nullptr)));
static_assert(noexcept(uno_create_model(nullptr, 0, nullptr, nullptr, 0)));
static_assert(noexcept(uno_create_unconstrained_model(nullptr, 0, 0)));
static_assert(noexcept(uno_set_constraints(nullptr, 0, nullptr, nullptr, nullptr, 0, nullptr, nullptr, nullptr)));
static_assert(noexcept(uno_set_lagrangian_hessian(nullptr, 0, 'L', nullptr, nullptr, nullptr)));
static_assert(noexcept(uno_create_solver()));
static_assert(noexcept(uno_set_solver_callbacks(nullptr, nullptr, nullptr, nullptr)));
static_assert(noexcept(uno_set_logger_stream_callback(nullptr, nullptr, nullptr)));
static_assert(noexcept(uno_reset_logger_stream(nullptr)));
static_assert(noexcept(uno_get_method_description(nullptr, nullptr)));
static_assert(noexcept(uno_get_optimization_status(nullptr, nullptr)));
static_assert(noexcept(uno_get_primal_solution_component(nullptr, 0, nullptr)));
static_assert(noexcept(uno_get_primal_solution(nullptr, nullptr)));
static_assert(noexcept(uno_get_number_iterations(nullptr, nullptr)));
static_assert(noexcept(uno_destroy_model(nullptr)));
static_assert(noexcept(uno_destroy_solver(nullptr)));

namespace {
   // owns a C API model of hs015 and a solver
   class CAPIResultTests: public ::testing::Test {
   protected:
      void SetUp() override {
         // a skip in SetUp skips the whole test
         if (available_preset() == nullptr) {
            GTEST_SKIP() << "no subproblem solver available";
         }
         this->model = create_hs015_model();
         this->solver = uno_create_solver();
         ASSERT_NE(this->model, nullptr);
         ASSERT_NE(this->solver, nullptr);
         ASSERT_EQ(uno_set_solver_string_option(this->solver, "logger", "SILENT"), UNO_OPTIONS_SUCCESS);
      }

      void TearDown() override {
         uno_destroy_solver(this->solver);
         uno_destroy_model(this->model);
      }

      // to be called as ASSERT_TRUE(this->optimize()): a failed assertion in a helper does not stop the test
      [[nodiscard]] bool optimize() const {
         return uno_set_solver_preset(this->solver, available_preset()) == UNO_OPTIONS_SUCCESS &&
            uno_optimize(this->solver, this->model) == UNO_OPTIMIZE_SUCCESS;
      }

      // checks that all the getters of the result return "expected_status"
      void expect_all_getters(uno_int expected_status) const {
         uno_int integer_value{-1};
         double double_value{-1.};
         std::vector<double> vector(2, -1.);
         EXPECT_EQ(uno_get_optimization_status(this->solver, &integer_value), expected_status);
         EXPECT_EQ(uno_get_solution_status(this->solver, &integer_value), expected_status);
         EXPECT_EQ(uno_get_solution_objective(this->solver, &double_value), expected_status);
         EXPECT_EQ(uno_get_primal_solution_component(this->solver, 0, &double_value), expected_status);
         EXPECT_EQ(uno_get_constraint_dual_solution_component(this->solver, 0, &double_value), expected_status);
         EXPECT_EQ(uno_get_lower_bound_dual_solution_component(this->solver, 0, &double_value), expected_status);
         EXPECT_EQ(uno_get_upper_bound_dual_solution_component(this->solver, 0, &double_value), expected_status);
         EXPECT_EQ(uno_get_primal_solution(this->solver, vector.data()), expected_status);
         EXPECT_EQ(uno_get_constraint_dual_solution(this->solver, vector.data()), expected_status);
         EXPECT_EQ(uno_get_lower_bound_dual_solution(this->solver, vector.data()), expected_status);
         EXPECT_EQ(uno_get_upper_bound_dual_solution(this->solver, vector.data()), expected_status);
         EXPECT_EQ(uno_get_solution_constraints(this->solver, vector.data()), expected_status);
         EXPECT_EQ(uno_get_solution_primal_feasibility(this->solver, &double_value), expected_status);
         EXPECT_EQ(uno_get_solution_stationarity(this->solver, &double_value), expected_status);
         EXPECT_EQ(uno_get_solution_complementarity(this->solver, &double_value), expected_status);
         EXPECT_EQ(uno_get_number_iterations(this->solver, &integer_value), expected_status);
         EXPECT_EQ(uno_get_cpu_time(this->solver, &double_value), expected_status);
         EXPECT_EQ(uno_get_number_objective_evaluations(this->solver, &integer_value), expected_status);
         EXPECT_EQ(uno_get_number_constraint_evaluations(this->solver, &integer_value), expected_status);
         EXPECT_EQ(uno_get_number_objective_gradient_evaluations(this->solver, &integer_value), expected_status);
         EXPECT_EQ(uno_get_number_jacobian_evaluations(this->solver, &integer_value), expected_status);
         EXPECT_EQ(uno_get_number_hessian_evaluations(this->solver, &integer_value), expected_status);
         EXPECT_EQ(uno_get_number_subproblem_solved_evaluations(this->solver, &integer_value), expected_status);
         const char* description{nullptr};
         EXPECT_EQ(uno_get_method_description(this->solver, &description), expected_status);
         if (expected_status != UNO_RESULT_SUCCESS) {
            // the values are untouched
            EXPECT_EQ(integer_value, -1);
            EXPECT_EQ(double_value, -1.);
            EXPECT_EQ(vector, std::vector<double>(2, -1.));
            EXPECT_EQ(description, nullptr);
         }
      }

      uno_model* model{nullptr};
      uno_solver* solver{nullptr};
   };
} // namespace

TEST_F(CAPIResultTests, ResultIsNotAvailableBeforeTheSolve) {
   this->expect_all_getters(UNO_RESULT_NOT_AVAILABLE);
}

TEST_F(CAPIResultTests, ResultIsAvailableAfterTheSolve) {
   ASSERT_TRUE(this->optimize());
   this->expect_all_getters(UNO_RESULT_SUCCESS);
}

TEST_F(CAPIResultTests, ResultIsDiscardedByAFailedOptimization) {
   ASSERT_TRUE(this->optimize());
   EXPECT_EQ(uno_optimize(this->solver, nullptr), UNO_OPTIMIZE_INVALID_ARGUMENT);
   this->expect_all_getters(UNO_RESULT_NOT_AVAILABLE);
}

TEST_F(CAPIResultTests, NullArgumentsAreInvalid) {
   ASSERT_TRUE(this->optimize());
   uno_int integer_value{-1};
   double double_value{-1.};
   EXPECT_EQ(uno_get_optimization_status(nullptr, &integer_value), UNO_RESULT_INVALID_ARGUMENT);
   EXPECT_EQ(uno_get_optimization_status(this->solver, nullptr), UNO_RESULT_INVALID_ARGUMENT);
   EXPECT_EQ(uno_get_solution_objective(nullptr, &double_value), UNO_RESULT_INVALID_ARGUMENT);
   EXPECT_EQ(uno_get_solution_objective(this->solver, nullptr), UNO_RESULT_INVALID_ARGUMENT);
   EXPECT_EQ(uno_get_primal_solution_component(nullptr, 0, &double_value), UNO_RESULT_INVALID_ARGUMENT);
   EXPECT_EQ(uno_get_primal_solution_component(this->solver, 0, nullptr), UNO_RESULT_INVALID_ARGUMENT);
   EXPECT_EQ(uno_get_primal_solution(nullptr, &double_value), UNO_RESULT_INVALID_ARGUMENT);
   EXPECT_EQ(uno_get_primal_solution(this->solver, nullptr), UNO_RESULT_INVALID_ARGUMENT);
   EXPECT_EQ(uno_get_number_iterations(this->solver, nullptr), UNO_RESULT_INVALID_ARGUMENT);
   EXPECT_EQ(integer_value, -1);
   EXPECT_EQ(double_value, -1.);
}

TEST_F(CAPIResultTests, SolutionValues) {
   ASSERT_TRUE(this->optimize());
   uno_int optimization_status{-1};
   uno_int solution_status{-1};
   double objective{};
   EXPECT_EQ(uno_get_optimization_status(this->solver, &optimization_status), UNO_RESULT_SUCCESS);
   EXPECT_EQ(optimization_status, UNO_SUCCESS);
   EXPECT_EQ(uno_get_solution_status(this->solver, &solution_status), UNO_RESULT_SUCCESS);
   EXPECT_EQ(solution_status, UNO_FEASIBLE_KKT_POINT);
   EXPECT_EQ(uno_get_solution_objective(this->solver, &objective), UNO_RESULT_SUCCESS);
   EXPECT_NEAR(objective, 306.5, 1e-6 * 306.5);

   // the vectors and their components agree
   std::vector<double> primal_solution(number_variables);
   ASSERT_EQ(uno_get_primal_solution(this->solver, primal_solution.data()), UNO_RESULT_SUCCESS);
   EXPECT_NEAR(primal_solution[0], 0.5, 1e-6);
   EXPECT_NEAR(primal_solution[1], 2., 1e-6);
   for (uno_int index = 0; index < number_variables; ++index) {
      double component{};
      EXPECT_EQ(uno_get_primal_solution_component(this->solver, index, &component), UNO_RESULT_SUCCESS);
      EXPECT_EQ(component, primal_solution[static_cast<size_t>(index)]);
   }
   std::vector<double> constraint_dual_solution(number_constraints);
   ASSERT_EQ(uno_get_constraint_dual_solution(this->solver, constraint_dual_solution.data()), UNO_RESULT_SUCCESS);
   for (uno_int index = 0; index < number_constraints; ++index) {
      double component{};
      EXPECT_EQ(uno_get_constraint_dual_solution_component(this->solver, index, &component), UNO_RESULT_SUCCESS);
      EXPECT_EQ(component, constraint_dual_solution[static_cast<size_t>(index)]);
   }

   uno_int number_iterations{};
   EXPECT_EQ(uno_get_number_iterations(this->solver, &number_iterations), UNO_RESULT_SUCCESS);
   EXPECT_GT(number_iterations, 0);
   uno_int number_evaluations{};
   EXPECT_EQ(uno_get_number_objective_evaluations(this->solver, &number_evaluations), UNO_RESULT_SUCCESS);
   EXPECT_GT(number_evaluations, 0);
   double cpu_time{-1.};
   EXPECT_EQ(uno_get_cpu_time(this->solver, &cpu_time), UNO_RESULT_SUCCESS);
   EXPECT_GE(cpu_time, 0.);
}

TEST_F(CAPIResultTests, InvalidIndicesAreRejected) {
   ASSERT_TRUE(this->optimize());
   double value{-1.};
   EXPECT_EQ(uno_get_primal_solution_component(this->solver, -1, &value), UNO_RESULT_INVALID_INDEX);
   EXPECT_EQ(uno_get_primal_solution_component(this->solver, number_variables, &value), UNO_RESULT_INVALID_INDEX);
   EXPECT_EQ(uno_get_constraint_dual_solution_component(this->solver, number_constraints, &value), UNO_RESULT_INVALID_INDEX);
   EXPECT_EQ(uno_get_lower_bound_dual_solution_component(this->solver, number_variables, &value), UNO_RESULT_INVALID_INDEX);
   EXPECT_EQ(uno_get_upper_bound_dual_solution_component(this->solver, -1, &value), UNO_RESULT_INVALID_INDEX);
   EXPECT_EQ(value, -1.);
}

TEST_F(CAPIResultTests, MethodDescription) {
   const char* description{nullptr};
   // like the other getters, no description before the solve
   EXPECT_EQ(uno_get_method_description(this->solver, &description), UNO_RESULT_NOT_AVAILABLE);
   EXPECT_EQ(description, nullptr);
   EXPECT_EQ(uno_get_method_description(nullptr, &description), UNO_RESULT_INVALID_ARGUMENT);
   EXPECT_EQ(uno_get_method_description(this->solver, nullptr), UNO_RESULT_INVALID_ARGUMENT);
   ASSERT_TRUE(this->optimize());
   ASSERT_EQ(uno_get_method_description(this->solver, &description), UNO_RESULT_SUCCESS);
   const std::string method{description};
   EXPECT_TRUE(method.find("SQP") != std::string::npos || method.find("interior-point") != std::string::npos) << method;
}

// an unavailable result is not an error: no message is printed, so that the getters can be used to check whether a
// result is available. The messages printed outside of an optimization go to the global logger stream
TEST_F(CAPIResultTests, NotAvailableResultPrintsNothing) {
   const Level level = Logger::level;
   std::ostream* const stream = Logger::stream;
   Logger::set_logger("WARNING");
   std::ostringstream messages{};
   Logger::set_stream(messages);
   double value{-1.};
   EXPECT_EQ(uno_get_solution_objective(this->solver, &value), UNO_RESULT_NOT_AVAILABLE);
   const std::string messages_not_available = messages.str();
   // an invalid argument is reported
   EXPECT_EQ(uno_get_solution_objective(nullptr, &value), UNO_RESULT_INVALID_ARGUMENT);
   const std::string messages_invalid_argument = messages.str();
   Logger::stream = stream;
   Logger::level = level;
   EXPECT_EQ(messages_not_available, "");
   EXPECT_NE(messages_invalid_argument, "");
}

TEST(CAPIModel, InvalidArgumentsAreRejected) {
   uno_int major{-1};
   uno_get_version(&major, nullptr, nullptr);
   EXPECT_EQ(major, UNO_VERSION_MAJOR);
   EXPECT_EQ(uno_create_unconstrained_model(nullptr, 2, UNO_ZERO_BASED_INDEXING), nullptr);
   EXPECT_EQ(uno_create_unconstrained_model("MINLP", 2, UNO_ZERO_BASED_INDEXING), nullptr);
   EXPECT_EQ(uno_create_unconstrained_model(UNO_PROBLEM_NONLINEAR, 0, UNO_ZERO_BASED_INDEXING), nullptr);
   EXPECT_EQ(uno_create_model(nullptr, 2, nullptr, nullptr, UNO_ZERO_BASED_INDEXING), nullptr);
   EXPECT_EQ(uno_create_model(UNO_PROBLEM_NONLINEAR, 2, nullptr, nullptr, 2), nullptr);

   uno_model* model = uno_create_unconstrained_model(UNO_PROBLEM_NONLINEAR, 2, UNO_ZERO_BASED_INDEXING);
   ASSERT_NE(model, nullptr);
   // null sparsity patterns
   EXPECT_NE(uno_set_constraints(model, 2, constraint_functions, nullptr, nullptr, 4, nullptr, nullptr, jacobian), UNO_MODEL_SUCCESS);
   EXPECT_NE(uno_set_constraints(model, 2, constraint_functions, nullptr, nullptr, -1, nullptr, nullptr, jacobian), UNO_MODEL_SUCCESS);
   EXPECT_NE(uno_set_lagrangian_hessian(model, 3, UNO_LOWER_TRIANGLE, nullptr, nullptr, lagrangian_hessian), UNO_MODEL_SUCCESS);
   EXPECT_NE(uno_set_variable_lower_bound(model, 2, 0.), UNO_MODEL_SUCCESS);
   EXPECT_NE(uno_set_objective(nullptr, UNO_MINIMIZE, objective_function, objective_gradient), UNO_MODEL_SUCCESS);
   EXPECT_EQ(uno_set_objective(model, UNO_MINIMIZE, objective_function, objective_gradient), UNO_MODEL_SUCCESS);
   uno_destroy_model(model);
}

TEST(CAPISolver, CallbacksAndLoggerStream) {
   uno_solver* solver = uno_create_solver();
   ASSERT_NE(solver, nullptr);
   EXPECT_NE(uno_set_solver_callbacks(nullptr, nullptr, nullptr, nullptr), UNO_SOLVER_SUCCESS);
   EXPECT_EQ(uno_set_solver_callbacks(solver, nullptr, nullptr, nullptr), UNO_SOLVER_SUCCESS);
   EXPECT_EQ(uno_set_logger_stream_callback(solver, nullptr, nullptr), UNO_SOLVER_SUCCESS);
   EXPECT_EQ(uno_reset_logger_stream(solver), UNO_SOLVER_SUCCESS);
   uno_destroy_solver(solver);
}
