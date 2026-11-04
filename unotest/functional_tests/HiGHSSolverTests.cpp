// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <chrono>
#include <limits>
#include <stdexcept>
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

// Bug: HiGHS' QP solver may return a primal solution that is inconsistent with its basis: below (a QP subproblem of the
// bound-constrained Rosenbrock problem near the solution), d0 is reported nonbasic at its upper bound 1.19e-5, with a
// nonzero bound multiplier, but its value is 0. Uno used the values as is: the SQP iterates never reached the active
// bound, the complementarity residual stalled at 1.2e-5 and Uno made null steps until the iteration limit.
TEST(HiGHSSolver, NonbasicVariablesLieOnTheirActiveBound) {
   Options options;
   options.set_bool("print_subproblem", false);
   HiGHSSolver solver(options);
   const double upper_bound = 1.19048176042136e-05;
   const Vector<double> linear_objective{-1.00002380963521, 5.55111512312578e-15};
   const Vector<uno_int> jacobian_row_indices{};
   const Vector<uno_int> jacobian_column_indices{};
   const Vector<double> jacobian_values{};
   const Vector<uno_int> hessian_row_indices{0, 1, 1};
   const Vector<uno_int> hessian_column_indices{0, 0, 1};
   const Vector<double> hessian_values{201.990476259296, -199.9952380729585, 200.};
   solver.get_quadratic_program().fill(linear_objective, jacobian_row_indices, jacobian_column_indices, jacobian_values,
      hessian_row_indices, hessian_column_indices, hessian_values, {-0.316712890936696, -0.316712890936696},
      {upper_bound, 0.316712890936696}, {}, {});
   Direction direction(2, 0);
   const WarmstartInformation warmstart_information{};
   Statistics statistics;
   solver.solve(statistics, Vector<double>(2), direction, warmstart_information);
   ASSERT_EQ(direction.status, SubproblemStatus::OPTIMAL);

   // the upper bound of d0 is active (the unconstrained minimizer is d0 = 0.5): complementarity must hold
   EXPECT_NEAR(direction.multipliers.upper_bounds[0], -1., 1e-3);
   EXPECT_DOUBLE_EQ(direction.primals[0], upper_bound);
}
namespace {
   // min 1/2 (x0^2 + x1^2) - x0 x1 / 2 - x0  s.t.  x0 + x1 <= 1,  0 <= x <= 10
   struct SmallQP {
      Vector<double> linear_objective{-1., 0.};
      Vector<uno_int> jacobian_row_indices{0, 0};
      Vector<uno_int> jacobian_column_indices{0, 1};
      Vector<double> jacobian_values{1., 1.};
      Vector<uno_int> hessian_row_indices{0, 1, 1};
      Vector<uno_int> hessian_column_indices{0, 0, 1};
      Vector<double> hessian_values{1., -0.5, 1.};
      std::vector<double> variables_lower_bounds{0., 0.};
      std::vector<double> variables_upper_bounds{10., 10.};
      std::vector<double> constraints_lower_bounds{-Inf};
      std::vector<double> constraints_upper_bounds{1.};

      void fill(QuadraticProgram& quadratic_program) const {
         quadratic_program.fill(linear_objective, jacobian_row_indices, jacobian_column_indices, jacobian_values,
            hessian_row_indices, hessian_column_indices, hessian_values, variables_lower_bounds, variables_upper_bounds,
            constraints_lower_bounds, constraints_upper_bounds);
      }
   };

   Direction solve(HiGHSSolver& solver) {
      Direction direction(2, 1);
      const WarmstartInformation warmstart_information{};
      Statistics statistics;
      solver.solve(statistics, Vector<double>(2), direction, warmstart_information);
      return direction;
   }
} // namespace

// Bug: the COO -> CSC conversion trusted the indices: a column index >= number of variables advanced the column
// counter past the end of the column starts array (out-of-bounds write), and inconsistent array lengths were read out
// of bounds. Malformed data must be rejected before anything is written.
TEST(HiGHSSolver, MalformedCOODataIsRejected) {
   Options options;
   options.set_bool("print_subproblem", false);
   HiGHSSolver solver(options);
   QuadraticProgram& quadratic_program = solver.get_quadratic_program();
   {
      SmallQP qp;
      qp.jacobian_column_indices[1] = 2; // out of range
      EXPECT_THROW(qp.fill(quadratic_program), std::invalid_argument);
   }
   {
      SmallQP qp;
      qp.jacobian_row_indices[0] = -1; // out of range
      EXPECT_THROW(qp.fill(quadratic_program), std::invalid_argument);
   }
   {
      SmallQP qp;
      qp.hessian_row_indices[2] = 5; // out of range
      EXPECT_THROW(qp.fill(quadratic_program), std::invalid_argument);
   }
   {
      SmallQP qp;
      qp.jacobian_values = Vector<double>{1.}; // inconsistent lengths
      EXPECT_THROW(qp.fill(quadratic_program), std::invalid_argument);
   }
   {
      SmallQP qp;
      qp.hessian_values[0] = std::numeric_limits<double>::quiet_NaN();
      EXPECT_THROW(qp.fill(quadratic_program), std::invalid_argument);
   }
   {
      SmallQP qp;
      qp.variables_upper_bounds.pop_back(); // inconsistent bound size
      EXPECT_THROW(qp.fill(quadratic_program), std::invalid_argument);
   }
   // valid data can still be filled and solved afterwards
   const SmallQP qp;
   qp.fill(quadratic_program);
   EXPECT_EQ(solve(solver).status, SubproblemStatus::OPTIMAL);
}

// the COO Hessian may be given as either triangle: HiGHS expects the lower one
TEST(HiGHSSolver, UpperTriangularHessianIsMirrored) {
   Options options;
   options.set_bool("print_subproblem", false);
   HiGHSSolver lower_solver(options);
   HiGHSSolver upper_solver(options);
   const SmallQP lower_qp;
   SmallQP upper_qp;
   upper_qp.hessian_row_indices = Vector<uno_int>{0, 0, 1};
   upper_qp.hessian_column_indices = Vector<uno_int>{0, 1, 1};
   lower_qp.fill(lower_solver.get_quadratic_program());
   upper_qp.fill(upper_solver.get_quadratic_program());
   const Direction lower_direction = solve(lower_solver);
   const Direction upper_direction = solve(upper_solver);
   ASSERT_EQ(lower_direction.status, SubproblemStatus::OPTIMAL);
   ASSERT_EQ(upper_direction.status, SubproblemStatus::OPTIMAL);
   for (size_t index: Range(2)) {
      EXPECT_NEAR(lower_direction.primals[index], upper_direction.primals[index], 1e-8);
   }
}
