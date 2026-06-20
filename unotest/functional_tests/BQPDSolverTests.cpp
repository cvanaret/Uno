// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "ingredients/subproblem_solvers/BQPD/BQPDSolver.hpp"
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
// BQPDQuadraticProgram::build() converts the COO matrices into BQPD's native (weak-CSR) layout internally.

TEST(BQPDSolver, LP) {
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
   options.set_string("BQPD_kmax_heuristic", "filtersqp");
   BQPDSolver solver(options);

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
   const std::vector<double> variables_upper_bounds{4., INF<double>};
   const std::vector<double> constraints_lower_bounds{-INF<double>, 5., 6.};
   const std::vector<double> constraints_upper_bounds{7., 15., INF<double>};

   QuadraticProgram& quadratic_program = solver.get_quadratic_program();
   quadratic_program.build(linear_objective, jacobian_row_indices, jacobian_column_indices, jacobian_values,
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

TEST(BQPDSolver, QP) {
   // https://doc.cgal.org/latest/QP_solver/index.html#title4
   // Min    f  =  1/2 * (2 x_0^2 + 8 x_1^2) - 32 x_1
   // s.t.         x_0 + x_1 <= 7
   //             -x_0 + 2x_1 <= 4
   // 0 <= x_0; 0 <= x_1 <= 4
   const size_t number_variables = 2;
   const size_t number_constraints = 2;

   Options options;
   options.set_bool("print_subproblem", false);
   options.set_string("BQPD_kmax_heuristic", "filtersqp");
   BQPDSolver solver(options);

   // dense objective gradient (linear part)
   const Vector<double> linear_objective{0., -32.};
   // COO constraint Jacobian
   const Vector<uno_int> jacobian_row_indices{0, 0, 1, 1};
   const Vector<uno_int> jacobian_column_indices{0, 1, 0, 1};
   const Vector<double> jacobian_values{1., 1., -1., 2.};
   // COO Lagrangian Hessian diag(2, 8) (one triangle)
   const Vector<uno_int> hessian_row_indices{0, 1};
   const Vector<uno_int> hessian_column_indices{0, 1};
   const Vector<double> hessian_values{2., 8.};
   // bounds
   const std::vector<double> variables_lower_bounds{0., 0.};
   const std::vector<double> variables_upper_bounds{INF<double>, 4.};
   const std::vector<double> constraints_lower_bounds{-INF<double>, -INF<double>};
   const std::vector<double> constraints_upper_bounds{7., 4.};

   QuadraticProgram& quadratic_program = solver.get_quadratic_program();
   quadratic_program.build(linear_objective, jacobian_row_indices, jacobian_column_indices, jacobian_values,
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
   const std::vector<double> primals_reference{2., 3.};
   for (size_t index: Range(number_variables)) {
      EXPECT_NEAR(direction.primals[index], primals_reference[index], tolerance);
   }
   // check duals
   const std::vector<double> constraint_duals_reference{0., -4.};
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
