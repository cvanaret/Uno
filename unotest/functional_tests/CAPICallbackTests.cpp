// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Tests of the user callbacks of the C API: contract of the return codes, and what the callbacks receive

#include <gtest/gtest.h>
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
