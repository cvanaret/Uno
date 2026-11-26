// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Signatures of the evaluation callbacks of the C API:
// - the callbacks that depend on the constraints all receive number_variables and number_constraints first;
// - every callback receives a flag "new_x" right after "x" (as in IPOPT): true if x differs from the point of the
// previous evaluation (of any function), false otherwise. The values that depend only on x can then be reused;
// - the operator callbacks have no "evaluate_at_x" flag (Uno always passed true)

#include <gtest/gtest.h>
#include <cmath>
#include <type_traits>
#include <vector>
#include "CAPITestProblem.hpp"

using namespace uno;
using namespace uno::capi_test;

// (n, x, new_x, f, user_data)
static_assert(std::is_same_v<uno_objective_callback, uno_int (*)(uno_int, const double*, bool, double*, void*)>);
// (n, x, new_x, gradient, user_data)
static_assert(std::is_same_v<uno_objective_gradient_callback, uno_int (*)(uno_int, const double*, bool, double*, void*)>);
// (n, m, x, new_x, c, user_data)
static_assert(std::is_same_v<uno_constraints_callback,
   uno_int (*)(uno_int, uno_int, const double*, bool, double*, void*)>);
// (n, m, nnz, x, new_x, values, user_data)
static_assert(std::is_same_v<uno_constraints_jacobian_callback,
   uno_int (*)(uno_int, uno_int, uno_int, const double*, bool, double*, void*)>);
// (n, m, nnz, x, new_x, rho, y, values, user_data)
static_assert(std::is_same_v<uno_lagrangian_hessian_callback,
   uno_int (*)(uno_int, uno_int, uno_int, const double*, bool, double, const double*, double*, void*)>);
// (n, m, x, new_x, vector, result, user_data)
static_assert(std::is_same_v<uno_constraints_jacobian_operator_callback,
   uno_int (*)(uno_int, uno_int, const double*, bool, const double*, double*, void*)>);
static_assert(std::is_same_v<uno_constraints_jacobian_transposed_operator_callback,
   uno_int (*)(uno_int, uno_int, const double*, bool, const double*, double*, void*)>);
// (n, m, x, new_x, rho, y, vector, result, user_data)
static_assert(std::is_same_v<uno_lagrangian_hessian_operator_callback,
   uno_int (*)(uno_int, uno_int, const double*, bool, double, const double*, const double*, double*, void*)>);

namespace {
   // every evaluation: the point and the flag
   struct Evaluation {
      std::vector<double> x;
      bool new_x;
   };

   struct EvaluationLog {
      std::vector<Evaluation> evaluations;

      void record(const double* x, bool new_x) {
         this->evaluations.push_back({{x, x + number_variables}, new_x});
      }
   };

   EvaluationLog& log_of(void* user_data) {
      return *static_cast<EvaluationLog*>(user_data);
   }

   uno_int logged_objective(uno_int n, const double* x, bool new_x, double* objective_value, void* user_data) {
      log_of(user_data).record(x, new_x);
      return objective_function(n, x, new_x, objective_value, user_data);
   }

   uno_int logged_gradient(uno_int n, const double* x, bool new_x, double* gradient, void* user_data) {
      log_of(user_data).record(x, new_x);
      return objective_gradient(n, x, new_x, gradient, user_data);
   }

   uno_int logged_constraints(uno_int n, uno_int m, const double* x, bool new_x, double* values, void* user_data) {
      log_of(user_data).record(x, new_x);
      return constraint_functions(n, m, x, new_x, values, user_data);
   }

   uno_int logged_jacobian(uno_int n, uno_int m, uno_int nnz, const double* x, bool new_x, double* values, void* user_data) {
      log_of(user_data).record(x, new_x);
      return jacobian(n, m, nnz, x, new_x, values, user_data);
   }

   uno_int logged_hessian(uno_int n, uno_int m, uno_int nnz, const double* x, bool new_x, double objective_multiplier,
         const double* multipliers, double* values, void* user_data) {
      log_of(user_data).record(x, new_x);
      return lagrangian_hessian(n, m, nnz, x, new_x, objective_multiplier, multipliers, values, user_data);
   }

   uno_model* create_logged_hs015_model(EvaluationLog& log) {
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
      uno_set_objective(model, UNO_MINIMIZE, logged_objective, logged_gradient);
      uno_set_constraints(model, number_constraints, logged_constraints, constraints_lower_bounds,
         constraints_upper_bounds, 4, jacobian_row_indices, jacobian_column_indices, logged_jacobian);
      uno_set_lagrangian_hessian(model, 3, UNO_LOWER_TRIANGLE, hessian_row_indices, hessian_column_indices,
         logged_hessian);
      uno_set_initial_primal_iterate(model, x0);
      uno_set_user_data(model, &log);
      return model;
   }

   // new_x is true exactly when x differs from the point of the previous evaluation (true at the first evaluation)
   void expect_consistent_flags(const std::vector<Evaluation>& evaluations) {
      ASSERT_FALSE(evaluations.empty());
      EXPECT_TRUE(evaluations.front().new_x) << "the first evaluation of a solve is at a new point";
      size_t number_new{0}, number_same{0};
      for (size_t index = 1; index < evaluations.size(); ++index) {
         const bool same_point = (evaluations[index].x == evaluations[index - 1].x);
         EXPECT_EQ(evaluations[index].new_x, !same_point) << "at evaluation " << index;
         (same_point ? number_same : number_new)++;
      }
      // the solve evaluates several functions at the same point, and visits several points
      EXPECT_GT(number_same, 0);
      EXPECT_GT(number_new, 0);
   }

   class CAPINewXTests: public ::testing::Test {
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

      uno_solver* solver{nullptr};
   };
} // namespace

TEST_F(CAPINewXTests, NewXFlagsThePointChanges) {
   EvaluationLog log;
   uno_model* model = create_logged_hs015_model(log);
   ASSERT_EQ(uno_optimize(this->solver, model), UNO_OPTIMIZE_SUCCESS);
   uno_destroy_model(model);
   expect_consistent_flags(log.evaluations);
}

// a new solve starts at a new point, even from the same starting point: the user cache may be stale
TEST_F(CAPINewXTests, NewXIsTrueAtTheFirstEvaluationOfEachSolve) {
   EvaluationLog log;
   uno_model* model = create_logged_hs015_model(log);
   ASSERT_EQ(uno_optimize(this->solver, model), UNO_OPTIMIZE_SUCCESS);
   const size_t number_evaluations_first_solve = log.evaluations.size();
   ASSERT_EQ(uno_optimize(this->solver, model), UNO_OPTIMIZE_SUCCESS);
   uno_destroy_model(model);
   ASSERT_GT(log.evaluations.size(), number_evaluations_first_solve);
   EXPECT_TRUE(log.evaluations[number_evaluations_first_solve].new_x);
}
