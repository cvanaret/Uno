// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Compile-time tests of the C API handles: the model and the solver are distinct opaque types (swapped arguments do
// not compile), and the getters accept a const solver

#include <gtest/gtest.h>
#include <type_traits>
#include "../../interfaces/C/Uno_C_API.h"

// the constructors return typed handles
static_assert(std::is_same_v<decltype(uno_create_model(nullptr, 0, nullptr, nullptr, 0)), uno_model*>);
static_assert(std::is_same_v<decltype(uno_create_unconstrained_model(nullptr, 0, 0)), uno_model*>);
static_assert(std::is_same_v<decltype(uno_create_solver()), uno_solver*>);

// the model and the solver cannot be swapped
static_assert(std::is_invocable_v<decltype(&uno_optimize), uno_solver*, uno_model*>);
static_assert(!std::is_invocable_v<decltype(&uno_optimize), uno_model*, uno_solver*>);
static_assert(!std::is_invocable_v<decltype(&uno_set_objective), uno_solver*, uno_int, uno_objective_callback,
   uno_objective_gradient_callback>);
static_assert(!std::is_invocable_v<decltype(&uno_set_solver_preset), uno_model*, const char*>);
static_assert(!std::is_invocable_v<decltype(&uno_destroy_model), uno_solver*>);
static_assert(!std::is_invocable_v<decltype(&uno_destroy_solver), uno_model*>);

// the getters do not modify the solver
static_assert(std::is_invocable_v<decltype(&uno_get_optimization_status), const uno_solver*, uno_int*>);
static_assert(std::is_invocable_v<decltype(&uno_get_primal_solution), const uno_solver*, double*>);
static_assert(std::is_invocable_v<decltype(&uno_get_number_iterations), const uno_solver*, uno_int*>);
static_assert(std::is_invocable_v<decltype(&uno_get_method_description), const uno_solver*, const char**>);
static_assert(std::is_invocable_v<decltype(&uno_get_solver_option_type), const uno_solver*, const char*>);
static_assert(std::is_invocable_v<decltype(&uno_get_solver_double_option), const uno_solver*, const char*, double*>);
// the setters do
static_assert(!std::is_invocable_v<decltype(&uno_set_solver_double_option), const uno_solver*, const char*, double>);

TEST(CAPIHandles, CreateAndDestroy) {
   uno_model* model = uno_create_unconstrained_model(UNO_PROBLEM_NONLINEAR, 1, UNO_ZERO_BASED_INDEXING);
   uno_solver* solver = uno_create_solver();
   EXPECT_NE(model, nullptr);
   EXPECT_NE(solver, nullptr);
   uno_destroy_solver(solver);
   uno_destroy_model(model);
}
