// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// The model and solver setters of the C API return status codes (0 upon success), like the options, result and
// optimize functions

#include <gtest/gtest.h>
#include <type_traits>
#include "CAPITestProblem.hpp"

using namespace uno::capi_test;

#define UNO_RETURNS_STATUS(call) static_assert(std::is_same_v<decltype(call), uno_int>, #call " must return uno_int")

UNO_RETURNS_STATUS(uno_set_model_name(nullptr, nullptr));
UNO_RETURNS_STATUS(uno_set_variables_lower_bounds(nullptr, nullptr));
UNO_RETURNS_STATUS(uno_set_variables_upper_bounds(nullptr, nullptr));
UNO_RETURNS_STATUS(uno_set_variable_lower_bound(nullptr, 0, 0.));
UNO_RETURNS_STATUS(uno_set_variable_upper_bound(nullptr, 0, 0.));
UNO_RETURNS_STATUS(uno_set_objective(nullptr, 0, nullptr, nullptr));
UNO_RETURNS_STATUS(uno_set_constraints(nullptr, 0, nullptr, nullptr, nullptr, 0, nullptr, nullptr, nullptr));
UNO_RETURNS_STATUS(uno_set_constraints_lower_bounds(nullptr, nullptr));
UNO_RETURNS_STATUS(uno_set_constraints_upper_bounds(nullptr, nullptr));
UNO_RETURNS_STATUS(uno_set_constraint_lower_bound(nullptr, 0, 0.));
UNO_RETURNS_STATUS(uno_set_constraint_upper_bound(nullptr, 0, 0.));
UNO_RETURNS_STATUS(uno_set_jacobian_operator(nullptr, nullptr));
UNO_RETURNS_STATUS(uno_set_jacobian_transposed_operator(nullptr, nullptr));
UNO_RETURNS_STATUS(uno_set_lagrangian_hessian(nullptr, 0, 'L', nullptr, nullptr, nullptr));
UNO_RETURNS_STATUS(uno_set_lagrangian_hessian_operator(nullptr, nullptr));
UNO_RETURNS_STATUS(uno_set_lagrangian_sign_convention(nullptr, 0));
UNO_RETURNS_STATUS(uno_set_user_data(nullptr, nullptr));
UNO_RETURNS_STATUS(uno_set_initial_primal_iterate_component(nullptr, 0, 0.));
UNO_RETURNS_STATUS(uno_set_initial_dual_iterate_component(nullptr, 0, 0.));
UNO_RETURNS_STATUS(uno_set_initial_primal_iterate(nullptr, nullptr));
UNO_RETURNS_STATUS(uno_set_initial_dual_iterate(nullptr, nullptr));
UNO_RETURNS_STATUS(uno_set_solver_callbacks(nullptr, nullptr, nullptr, nullptr));

namespace {
   class CAPIStatusCodeTests: public ::testing::Test {
   protected:
      void SetUp() override {
         this->model = create_hs015_model();
         ASSERT_NE(this->model, nullptr);
      }

      void TearDown() override {
         uno_destroy_model(this->model);
      }

      uno_model* model{nullptr};
   };
} // namespace

TEST_F(CAPIStatusCodeTests, SuccessIsZero) {
   EXPECT_EQ(UNO_MODEL_SUCCESS, 0);
   const double bounds[] = {-1., -1.};
   const double x0[] = {0., 0.};
   EXPECT_EQ(uno_set_model_name(this->model, "hs015"), UNO_MODEL_SUCCESS);
   EXPECT_EQ(uno_set_variables_lower_bounds(this->model, bounds), UNO_MODEL_SUCCESS);
   EXPECT_EQ(uno_set_variable_upper_bound(this->model, 0, 0.5), UNO_MODEL_SUCCESS);
   EXPECT_EQ(uno_set_objective(this->model, UNO_MINIMIZE, objective_function, objective_gradient), UNO_MODEL_SUCCESS);
   EXPECT_EQ(uno_set_constraint_lower_bound(this->model, 1, 0.), UNO_MODEL_SUCCESS);
   EXPECT_EQ(uno_set_lagrangian_sign_convention(this->model, UNO_MULTIPLIER_POSITIVE), UNO_MODEL_SUCCESS);
   EXPECT_EQ(uno_set_user_data(this->model, nullptr), UNO_MODEL_SUCCESS);
   EXPECT_EQ(uno_set_initial_primal_iterate(this->model, x0), UNO_MODEL_SUCCESS);
   EXPECT_EQ(uno_set_initial_dual_iterate_component(this->model, 1, 0.), UNO_MODEL_SUCCESS);
}

TEST_F(CAPIStatusCodeTests, NullModelIsInvalidArgument) {
   EXPECT_EQ(uno_set_model_name(nullptr, "hs015"), UNO_MODEL_INVALID_ARGUMENT);
   EXPECT_EQ(uno_set_objective(nullptr, UNO_MINIMIZE, objective_function, objective_gradient), UNO_MODEL_INVALID_ARGUMENT);
   EXPECT_EQ(uno_set_variable_lower_bound(nullptr, 0, 0.), UNO_MODEL_INVALID_ARGUMENT);
   EXPECT_EQ(uno_set_user_data(nullptr, nullptr), UNO_MODEL_INVALID_ARGUMENT);
}

TEST_F(CAPIStatusCodeTests, NullArrayIsInvalidArgument) {
   EXPECT_EQ(uno_set_initial_primal_iterate(this->model, nullptr), UNO_MODEL_INVALID_ARGUMENT);
   EXPECT_EQ(uno_set_initial_dual_iterate(this->model, nullptr), UNO_MODEL_INVALID_ARGUMENT);
   EXPECT_EQ(uno_set_variables_lower_bounds(this->model, nullptr), UNO_MODEL_INVALID_ARGUMENT);
   EXPECT_EQ(uno_set_constraints_upper_bounds(this->model, nullptr), UNO_MODEL_INVALID_ARGUMENT);
   // a positive number of Jacobian nonzeros requires the sparsity pattern
   EXPECT_EQ(uno_set_constraints(this->model, 2, constraint_functions, nullptr, nullptr, 4, nullptr, nullptr, jacobian),
      UNO_MODEL_INVALID_ARGUMENT);
   EXPECT_EQ(uno_set_constraints(this->model, 0, constraint_functions, nullptr, nullptr, 0, nullptr, nullptr, jacobian),
      UNO_MODEL_INVALID_ARGUMENT);
}

TEST_F(CAPIStatusCodeTests, OutOfRangeIndexIsInvalidIndex) {
   EXPECT_EQ(uno_set_variable_lower_bound(this->model, 2, 0.), UNO_MODEL_INVALID_INDEX);
   EXPECT_EQ(uno_set_variable_upper_bound(this->model, -1, 0.), UNO_MODEL_INVALID_INDEX);
   EXPECT_EQ(uno_set_constraint_lower_bound(this->model, 2, 0.), UNO_MODEL_INVALID_INDEX);
   EXPECT_EQ(uno_set_initial_primal_iterate_component(this->model, 2, 0.), UNO_MODEL_INVALID_INDEX);
   EXPECT_EQ(uno_set_initial_dual_iterate_component(this->model, 2, 0.), UNO_MODEL_INVALID_INDEX);
}

TEST_F(CAPIStatusCodeTests, InvalidValueIsInvalidValue) {
   EXPECT_EQ(uno_set_objective(this->model, 2, objective_function, objective_gradient), UNO_MODEL_INVALID_VALUE);
   EXPECT_EQ(uno_set_lagrangian_sign_convention(this->model, 0), UNO_MODEL_INVALID_VALUE);
   uno_int indices[] = {0};
   EXPECT_EQ(uno_set_lagrangian_hessian(this->model, 1, 'X', indices, indices, lagrangian_hessian), UNO_MODEL_INVALID_VALUE);
}

TEST(CAPISolverStatusCodes, SolverCallbacks) {
   EXPECT_EQ(UNO_SOLVER_SUCCESS, 0);
   uno_solver* solver = uno_create_solver();
   ASSERT_NE(solver, nullptr);
   EXPECT_EQ(uno_set_solver_callbacks(nullptr, nullptr, nullptr, nullptr), UNO_SOLVER_INVALID_ARGUMENT);
   EXPECT_EQ(uno_set_solver_callbacks(solver, nullptr, nullptr, nullptr), UNO_SOLVER_SUCCESS);
   uno_destroy_solver(solver);
}
