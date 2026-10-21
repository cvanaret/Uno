// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <chrono>
#include <thread>
#include <gtest/gtest.h>
#if !defined(_WIN32)
#include <csignal>
#include <sys/wait.h>
#include <unistd.h>
#endif
#include "ingredients/subproblem_solvers/HiGHS/HiGHSSolver.hpp"
#include "ingredients/subproblem_solvers/QuadraticProgram.hpp"
#include "ingredients/subproblem_solvers/SubproblemStatus.hpp"
#include "../interfaces/C/uno_int.h"
#include "linear_algebra/Vector.hpp"
#include "optimization/Direction.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "symbolic/Range.hpp"
#include "tools/Infinity.hpp"
#include "tools/Statistics.hpp"

using namespace uno;

// The quadratic program is described directly from data (no Subproblem): a dense objective gradient, a COO
// constraint Jacobian (row = constraint, column = variable) and a COO Lagrangian Hessian (empty for an LP).
// HiGHSQuadraticProgram::build() converts the COO matrices into HiGHS' native CSC layout internally.

TEST(HiGHSSolver, LP) {
   // https://ergo-code.github.io/HiGHS/stable/interfaces/cpp/library/
   // Min    f  =  x_0 +  x_1 + 3
   // s.t.                x_1 <= 7
   //        5 <=  x_0 + 2x_1 <= 15
   //        6 <= 3x_0 + 2x_1
   // 0 <= x_0 <= 4; 1 <= x_1
   const size_t number_variables = 2;
   const size_t number_constraints = 3;

   Options options;
   options.set_bool("print_subproblem", false);
   HiGHSSolver solver(options);

   // dense objective gradient
   const Vector<double> linear_objective{1., 1.};
   // COO constraint Jacobian
   const Vector<uno_int> jacobian_row_indices{0, 1, 1, 2, 2};
   const Vector<uno_int> jacobian_column_indices{1, 0, 1, 0, 1};
   const Vector<double> jacobian_values{1., 1., 2., 3., 2.};
   // empty Hessian (LP)
   const Vector<uno_int> hessian_row_indices{};
   const Vector<uno_int> hessian_column_indices{};
   const Vector<double> hessian_values{};
   // bounds
   const std::vector<double> variables_lower_bounds{0., 1.};
   const std::vector<double> variables_upper_bounds{4., Inf};
   const std::vector<double> constraints_lower_bounds{-Inf, 5., 6.};
   const std::vector<double> constraints_upper_bounds{7., 15., Inf};

   QuadraticProgram& quadratic_program = solver.get_quadratic_program();
   quadratic_program.fill(linear_objective, jacobian_row_indices, jacobian_column_indices, jacobian_values,
      hessian_row_indices, hessian_column_indices, hessian_values,
      variables_lower_bounds, variables_upper_bounds, constraints_lower_bounds, constraints_upper_bounds);

   Direction direction(number_variables, number_constraints);
   WarmstartInformation warmstart_information{};
   const Vector<double> initial_point{0., 0.};
   Statistics statistics;
   solver.solve(statistics, initial_point, direction, warmstart_information);

   ASSERT_EQ(direction.status, SubproblemStatus::OPTIMAL);

   const double tolerance = 1e-8;
   // check primals
   const std::vector<double> primals_reference{0.5, 2.25};
   for (size_t index: Range(number_variables)) {
      EXPECT_NEAR(direction.primals[index], primals_reference[index], tolerance);
   }
   // check duals
   const std::vector<double> constraint_duals_reference{0., 0.25, 0.25};
   for (size_t index: Range(number_constraints)) {
      EXPECT_NEAR(direction.multipliers.constraints[index], constraint_duals_reference[index], tolerance);
   }
   const std::vector<double> lower_bound_duals_reference{0., 0.};
   const std::vector<double> upper_bound_duals_reference{0., 0.};
   for (size_t index: Range(number_variables)) {
      EXPECT_NEAR(direction.multipliers.lower_bounds[index], lower_bound_duals_reference[index], tolerance);
      EXPECT_NEAR(direction.multipliers.upper_bounds[index], upper_bound_duals_reference[index], tolerance);
   }
}

// Bug: HiGHS' active-set QP solver may stall on convex QPs with a semidefinite Hessian, such as the l1 feasibility
// subproblems (no curvature in the elastic variables). The QP below is the restoration subproblem of HS15 (trust region,
// primal regularization): HiGHS 1.15 performs millions of iterations without terminating. Uno set no iteration limit, so
// that the solve never returned. The QP is solved in a child process that is killed after a timeout.
TEST(HiGHSSolver, StallingQPSolveReturns) {
#if defined(_WIN32)
   GTEST_SKIP() << "requires fork()";
#else
   const pid_t pid = fork();
   ASSERT_NE(pid, -1);
   if (pid == 0) { // child process: solve the QP
      Options options;
      options.set_bool("print_subproblem", false);
      HiGHSSolver solver(options);
      // min p + n + 1e-4/2 (d0^2 + d1^2)
      // s.t. 0.5 d0 + p >= 0, d0 + d1 + n >= -0.25, -10 <= d0 <= 0.5, -10 <= d1 <= 10, p >= -1, n >= 0
      const Vector<double> linear_objective{0., 0., 1., 1.};
      const Vector<uno_int> jacobian_row_indices{0, 0, 1, 1, 1};
      const Vector<uno_int> jacobian_column_indices{0, 2, 0, 1, 3};
      const Vector<double> jacobian_values{0.5, 1., 1., 1., 1.};
      const Vector<uno_int> hessian_row_indices{0, 1};
      const Vector<uno_int> hessian_column_indices{0, 1};
      const Vector<double> hessian_values{1e-4, 1e-4};
      solver.get_quadratic_program().fill(linear_objective, jacobian_row_indices, jacobian_column_indices, jacobian_values,
         hessian_row_indices, hessian_column_indices, hessian_values, {-10., -10., -1., 0.}, {0.5, 10., Inf, Inf},
         {0., -0.25}, {Inf, Inf});
      Direction direction(4, 2);
      const WarmstartInformation warmstart_information{};
      Statistics statistics;
      solver.solve(statistics, Vector<double>(4), direction, warmstart_information);
      _exit(0);
   }
   // parent process: wait for the child at most 60 seconds
   int status = 0;
   bool returned = false;
   for (size_t attempt = 0; attempt < 600 && !returned; ++attempt) {
      returned = (waitpid(pid, &status, WNOHANG) == pid);
      if (!returned) {
         std::this_thread::sleep_for(std::chrono::milliseconds(100));
      }
   }
   if (!returned) {
      kill(pid, SIGKILL);
      waitpid(pid, &status, 0);
   }
   ASSERT_TRUE(returned) << "the HiGHS QP solve did not return";
   EXPECT_TRUE(WIFEXITED(status) && WEXITSTATUS(status) == 0);
#endif
}