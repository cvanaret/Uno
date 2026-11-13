// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Integration tests for Ctrl+C handling: SIGINT is raised from inside a user evaluation (i.e. deep in the inner loop)
// and must be honored at the end of the current outer iteration, returning the current iterate.

#include <gtest/gtest.h>
#include <cmath>
#include <csignal>
#include <vector>
#include "../../interfaces/C/Uno_C_API.h"
#include "tools/Infinity.hpp"

using namespace uno;

namespace {
   struct InterruptData {
      uno_int objective_evaluations{0};
      uno_int raise_at_evaluation{0}; // 0: never raise
      bool raised{false};
   };

   // problem hs015, unlikely to converge in a single outer iteration from x0 = (-2, 1)
   constexpr uno_int number_variables = 2;
   constexpr uno_int number_constraints = 2;

   uno_int objective_function(uno_int /*number_variables*/, const double* x, double* objective_value, void* user_data) {
      auto* data = static_cast<InterruptData*>(user_data);
      ++data->objective_evaluations;
      // raise only once: the handler reinstalls SIG_DFL, so a second raise would terminate the test process
      if (!data->raised && data->objective_evaluations == data->raise_at_evaluation) {
         data->raised = true;
         std::raise(SIGINT); // synchronous on all platforms
      }
      *objective_value = 100. * std::pow(x[1] - x[0]*x[0], 2.) + std::pow(1. - x[0], 2.);
      return 0;
   }

   uno_int objective_gradient(uno_int /*number_variables*/, const double* x, double* gradient, void* /*user_data*/) {
      gradient[0] = 400.*std::pow(x[0], 3.) - 400.*x[0]*x[1] + 2.*x[0] - 2.;
      gradient[1] = 200.*(x[1] - x[0]*x[0]);
      return 0;
   }

   uno_int constraint_functions(uno_int /*number_variables*/, uno_int /*number_constraints*/, const double* x,
         double* constraint_values, void* /*user_data*/) {
      constraint_values[0] = x[0]*x[1];
      constraint_values[1] = x[0] + x[1]*x[1];
      return 0;
   }

   uno_int jacobian(uno_int /*number_variables*/, uno_int /*number_jacobian_nonzeros*/, const double* x,
         double* jacobian_values, void* /*user_data*/) {
      jacobian_values[0] = x[1];
      jacobian_values[1] = 1.;
      jacobian_values[2] = x[0];
      jacobian_values[3] = 2.*x[1];
      return 0;
   }

   uno_int lagrangian_hessian(uno_int /*number_variables*/, uno_int /*number_constraints*/, uno_int /*number_hessian_nonzeros*/,
         const double* x, double objective_multiplier, const double* multipliers, double* hessian_values, void* /*user_data*/) {
      hessian_values[0] = objective_multiplier*(1200.*x[0]*x[0] - 400.*x[1] + 2.);
      hessian_values[1] = -400.*objective_multiplier*x[0] - multipliers[0];
      hessian_values[2] = 200.*objective_multiplier - 2.*multipliers[1];
      return 0;
   }

   // pick a preset whose subproblem solver was compiled in; nullptr if none
   const char* available_preset() {
#if defined(HAS_BQPD)
      return "filtersqp";
#elif defined(HAS_HSL) || defined(HAS_MA57) || defined(HAS_MA27) || defined(HAS_MUMPS) || defined(HAS_SPRAL)
      return "ipopt";
#elif defined(HAS_HIGHS)
      return "filterslp";
#else
      return nullptr;
#endif
   }

   // owns the C API model and solver
   class HS015 {
   public:
      explicit HS015(InterruptData& data) {
         double variables_lower_bounds[] = {-Inf, -Inf};
         double variables_upper_bounds[] = {0.5, Inf};
         double constraints_lower_bounds[] = {1., 0.};
         double constraints_upper_bounds[] = {Inf, Inf};
         uno_int jacobian_row_indices[] = {0, 1, 0, 1};
         uno_int jacobian_column_indices[] = {0, 0, 1, 1};
         uno_int hessian_row_indices[] = {0, 1, 1};
         uno_int hessian_column_indices[] = {0, 0, 1};
         double x0[] = {-2., 1.};

         this->model = uno_create_model(UNO_PROBLEM_NONLINEAR, number_variables, variables_lower_bounds,
            variables_upper_bounds, UNO_ZERO_BASED_INDEXING);
         uno_set_objective(this->model, UNO_MINIMIZE, objective_function, objective_gradient);
         uno_set_constraints(this->model, number_constraints, constraint_functions, constraints_lower_bounds,
            constraints_upper_bounds, 4, jacobian_row_indices, jacobian_column_indices, jacobian);
         uno_set_lagrangian_hessian(this->model, 3, UNO_LOWER_TRIANGLE, hessian_row_indices, hessian_column_indices,
            lagrangian_hessian);
         uno_set_lagrangian_sign_convention(this->model, UNO_MULTIPLIER_NEGATIVE);
         uno_set_initial_primal_iterate(this->model, x0);
         uno_set_user_data(this->model, &data);

         this->solver = uno_create_solver();
         uno_set_solver_preset(this->solver, available_preset());
         uno_set_solver_string_option(this->solver, "logger", "SILENT");
      }

      ~HS015() {
         uno_destroy_solver(this->solver);
         uno_destroy_model(this->model);
      }

      HS015(const HS015&) = delete;
      HS015& operator=(const HS015&) = delete;

      void optimize() {
         uno_optimize(this->solver, this->model);
      }

      void* model{nullptr};
      void* solver{nullptr};
   };

   volatile std::sig_atomic_t sentinel_calls = 0;
   extern "C" void interrupt_tests_sentinel_handler(int /*signal*/) {
      sentinel_calls = sentinel_calls + 1;
   }

   class InterruptTests: public ::testing::Test {
   protected:
      void SetUp() override {
         if (available_preset() == nullptr) {
            GTEST_SKIP() << "no subproblem solver available";
         }
      }
   };
} // namespace

// SIGINT raised in the inner loop of the first outer iteration terminates after that iteration with the current iterate
TEST_F(InterruptTests, TerminatesAtEndOfOuterIteration) {
   InterruptData data{};
   data.raise_at_evaluation = 2; // first evaluation after the initial point
   HS015 problem(data);
   problem.optimize();

   ASSERT_TRUE(data.raised);
   EXPECT_EQ(uno_get_optimization_status(problem.solver), UNO_USER_INTERRUPTION);
   EXPECT_EQ(uno_get_number_iterations(problem.solver), 1);

   std::vector<double> primals(number_variables);
   uno_get_primal_solution(problem.solver, primals.data());
   for (double primal: primals) {
      EXPECT_TRUE(std::isfinite(primal));
   }
   EXPECT_TRUE(std::isfinite(uno_get_solution_objective(problem.solver)));
}

// the inner loop is not aborted: evaluations continue after the signal until the outer loop polls the flag
TEST_F(InterruptTests, InnerLoopRunsToCompletion) {
   InterruptData data{};
   data.raise_at_evaluation = 1; // during initialization
   HS015 problem(data);
   problem.optimize();

   ASSERT_TRUE(data.raised);
   EXPECT_EQ(uno_get_optimization_status(problem.solver), UNO_USER_INTERRUPTION);
   EXPECT_EQ(uno_get_number_iterations(problem.solver), 1);
   EXPECT_GT(data.objective_evaluations, 1);
}

// an interrupted solve does not leak its flag into the next solve
TEST_F(InterruptTests, FlagIsResetBetweenSolves) {
   InterruptData data{};
   data.raise_at_evaluation = 2;
   HS015 problem(data);
   problem.optimize();
   ASSERT_EQ(uno_get_optimization_status(problem.solver), UNO_USER_INTERRUPTION);

   // same solver and model, no signal this time
   data.raise_at_evaluation = 0;
   problem.optimize();
   EXPECT_NE(uno_get_optimization_status(problem.solver), UNO_USER_INTERRUPTION);
   EXPECT_GT(uno_get_number_iterations(problem.solver), 1);
}

// the host's SIGINT handler is restored after the solve, even after an interrupt
TEST_F(InterruptTests, RestoresHostHandler) {
   const auto default_handler = std::signal(SIGINT, interrupt_tests_sentinel_handler);
   ASSERT_NE(default_handler, SIG_ERR);
   sentinel_calls = 0;

   {
      InterruptData data{};
      data.raise_at_evaluation = 2;
      HS015 problem(data);
      problem.optimize();
      ASSERT_EQ(uno_get_optimization_status(problem.solver), UNO_USER_INTERRUPTION);
   }
   EXPECT_EQ(sentinel_calls, 0); // the solver's handler intercepted the signal

   std::raise(SIGINT); // must reach the host handler, not kill the process
   EXPECT_EQ(sentinel_calls, 1);

   std::signal(SIGINT, default_handler);
}