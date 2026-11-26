// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Tests of the construction of a model through the C API

#include <gtest/gtest.h>
#include <optional>
#include "CAPITestProblem.hpp"

using namespace uno;
using namespace uno::capi_test;

namespace {
   // min (x - 3)^2 s.t. lb <= x <= ub (as a general constraint c(x) = x)
   uno_int shifted_square(uno_int /*n*/, const double* x, double* objective_value, void* /*user_data*/) {
      *objective_value = (x[0] - 3.)*(x[0] - 3.);
      return 0;
   }

   uno_int shifted_square_gradient(uno_int /*n*/, const double* x, double* gradient, void* /*user_data*/) {
      gradient[0] = 2.*(x[0] - 3.);
      return 0;
   }

   uno_int identity_constraint(uno_int /*n*/, uno_int /*m*/, const double* x, double* values, void* /*user_data*/) {
      values[0] = x[0];
      return 0;
   }

   uno_int identity_jacobian(uno_int /*n*/, uno_int /*nnz*/, const double* /*x*/, double* values, void* /*user_data*/) {
      values[0] = 1.;
      return 0;
   }

   uno_int shifted_square_hessian(uno_int /*n*/, uno_int /*m*/, uno_int /*nnz*/, const double* /*x*/,
         double objective_multiplier, const double* /*multipliers*/, double* values, void* /*user_data*/) {
      values[0] = 2.*objective_multiplier;
      return 0;
   }

   // solves the problem and returns the primal solution (nullopt if the solve fails)
   std::optional<double> solve_with_constraint_bounds(const double* lower_bound, const double* upper_bound) {
      double x0[] = {0.};
      uno_int zero[] = {0};
      uno_model* model = uno_create_model(UNO_PROBLEM_NONLINEAR, 1, nullptr, nullptr, UNO_ZERO_BASED_INDEXING);
      uno_set_objective(model, UNO_MINIMIZE, shifted_square, shifted_square_gradient);
      uno_set_constraints(model, 1, identity_constraint, lower_bound, upper_bound, 1, zero, zero, identity_jacobian);
      uno_set_lagrangian_hessian(model, 1, UNO_LOWER_TRIANGLE, zero, zero, shifted_square_hessian);
      uno_set_initial_primal_iterate(model, x0);
      uno_solver* solver = uno_create_solver();
      uno_set_solver_string_option(solver, "logger", "SILENT");
      uno_set_solver_preset(solver, available_preset());
      std::optional<double> solution{};
      uno_int status{-1};
      double x{0.};
      if (uno_optimize(solver, model) == UNO_OPTIMIZE_SUCCESS &&
            uno_get_optimization_status(solver, &status) == UNO_RESULT_SUCCESS && status == UNO_SUCCESS &&
            uno_get_primal_solution_component(solver, 0, &x) == UNO_RESULT_SUCCESS) {
         solution = x;
      }
      uno_destroy_solver(solver);
      uno_destroy_model(model);
      return solution;
   }

   class CAPIModelTests: public ::testing::Test {
   protected:
      void SetUp() override {
         if (available_preset() == nullptr) {
            GTEST_SKIP() << "no subproblem solver available";
         }
      }
   };
} // namespace

// NULL constraint bounds are infinite, as NULL variable bounds: the constraint is inactive
TEST_F(CAPIModelTests, NullConstraintBoundsAreInfinite) {
   const std::optional<double> x = solve_with_constraint_bounds(nullptr, nullptr);
   ASSERT_TRUE(x.has_value());
   EXPECT_NEAR(*x, 3., 1e-6);
}

TEST_F(CAPIModelTests, NullConstraintLowerBoundIsMinusInfinity) {
   // the constraint x <= -1 is active (a lower bound of 0 would make the problem infeasible)
   const double upper_bound[] = {-1.};
   const std::optional<double> x = solve_with_constraint_bounds(nullptr, upper_bound);
   ASSERT_TRUE(x.has_value());
   EXPECT_NEAR(*x, -1., 1e-6);
}

TEST_F(CAPIModelTests, NullConstraintUpperBoundIsPlusInfinity) {
   const double lower_bound[] = {1.};
   const std::optional<double> x = solve_with_constraint_bounds(lower_bound, nullptr);
   ASSERT_TRUE(x.has_value());
   EXPECT_NEAR(*x, 3., 1e-6);
}

// regression: finite bounds are taken into account
TEST_F(CAPIModelTests, FiniteConstraintBoundsAreEnforced) {
   const double lower_bound[] = {-1.};
   const double upper_bound[] = {2.};
   const std::optional<double> x = solve_with_constraint_bounds(lower_bound, upper_bound);
   ASSERT_TRUE(x.has_value());
   EXPECT_NEAR(*x, 2., 1e-6);
}
