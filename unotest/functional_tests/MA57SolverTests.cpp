// Copyright (c) 2018-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <random>
#include <vector>
#include "ingredients/subproblem_solvers/MA57/MA57Solver.hpp"
#include "LinearSolverTests.hpp"

namespace uno {
   template <>
   constexpr size_t number_pivoting_tolerance_increases<MA57Solver> = 3;

   INSTANTIATE_TYPED_TEST_SUITE_P(MA57, LinearSolverTest, MA57Solver);

   // KKT matrix [H A^T; A 0] with a tiny diagonal in H: the pivots chosen by the analysis are rejected by the numerical
   // pivoting, and the delayed pivots need more workspace than the analysis predicts. The factorization is then suspended,
   // the workspace enlarged and the factorization resumed
   TEST(MA57Solver, WorkspaceEnlargement) {
      constexpr size_t number_variables = 300, number_constraints = 200;
      constexpr size_t dimension = number_variables + number_constraints;
      std::mt19937 generator(1);
      std::uniform_int_distribution<size_t> random_variable(0, number_variables - 1);
      std::uniform_real_distribution<double> random_value(-1., 1.);

      // lower triangle, Fortran-based indexing
      std::vector<uno_int> rows, columns;
      std::vector<double> values;
      const auto add_entry = [&](size_t row, size_t column, double value) {
         rows.push_back(static_cast<uno_int>(row + 1));
         columns.push_back(static_cast<uno_int>(column + 1));
         values.push_back(value);
      };
      for (size_t variable_index = 0; variable_index < number_variables; ++variable_index) {
         add_entry(variable_index, variable_index, 1e-14 * (1. + random_value(generator)));
      }
      for (size_t constraint_index = 0; constraint_index < number_constraints; ++constraint_index) {
         for (size_t term = 0; term < 6; ++term) {
            add_entry(number_variables + constraint_index, random_variable(generator), random_value(generator));
         }
      }

      MA57Solver solver;
      COOLinearSystem& linear_system = solver.get_coo_linear_system();
      linear_system.dimension = dimension;
      linear_system.number_nonzeros = values.size();
      linear_system.matrix_row_indices = rows;
      linear_system.matrix_column_indices = columns;
      linear_system.matrix_values = Vector<double>(values.size());
      for (size_t index = 0; index < values.size(); ++index) {
         linear_system.matrix_values[index] = values[index];
      }
      // right-hand side of the solution x = (1, ..., 1)
      linear_system.rhs = Vector<double>(dimension, 0.);
      for (size_t index = 0; index < values.size(); ++index) {
         const auto row = static_cast<size_t>(rows[index] - 1), column = static_cast<size_t>(columns[index] - 1);
         linear_system.rhs[row] += values[index];
         if (row != column) {
            linear_system.rhs[column] += values[index];
         }
      }
      solver.initialize_memory();
      solver.do_symbolic_analysis();
      solver.do_numerical_factorization(false);
      ASSERT_FALSE(solver.matrix_is_singular());

      Vector<double> solution(dimension);
      solution.fill(0.);
      solver.solve_indefinite_system(solution.data());
      // residual of the system
      std::vector<double> residual(linear_system.rhs.begin(), linear_system.rhs.end());
      for (size_t index = 0; index < values.size(); ++index) {
         const auto row = static_cast<size_t>(rows[index] - 1), column = static_cast<size_t>(columns[index] - 1);
         residual[row] -= values[index] * solution[column];
         if (row != column) {
            residual[column] -= values[index] * solution[row];
         }
      }
      double residual_norm = 0.;
      for (const double entry: residual) {
         residual_norm = std::max(residual_norm, std::abs(entry));
      }
      EXPECT_LT(residual_norm, 1e-8);
   }
}
