// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// The elimination of the fixed variables is transparent: the KKT point is returned in the space of the original model

#include <algorithm>
#include <string>
#include <vector>
#include <gtest/gtest.h>
#include "../ReductionTestProblems.hpp"
#include "Uno.hpp"
#include "ingredients/subproblem_solvers/QPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "optimization/Result.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"
#include "tools/Logger.hpp"

using namespace uno;

namespace {
   struct Solution {
      std::vector<double> primals, constraint_multipliers, lower_bound_multipliers, upper_bound_multipliers;
      double objective;
   };

   bool is_preset_available(const std::string& preset) {
      if (preset == "filtersqp") {
         return std::any_of(QPSolverFactory::available_solvers.begin(), QPSolverFactory::available_solvers.end(),
            [](const char* solver) { return std::string(solver) == "BQPD"; });
      }
      return !SymmetricIndefiniteLinearSolverFactory::available_solvers().empty();
   }

   Result solve(const Model& model, const std::string& preset, int max_iterations = -1) {
      Options options;
      DefaultOptions::load(options);
      Presets::set(options, preset);
      if (0 <= max_iterations) {
         options.set_integer("max_iterations", max_iterations);
      }
      Logger::level = SILENT;
      Uno uno;
      return uno.solve(model, options);
   }

   void expect_near(const Vector<double>& actual, const std::vector<double>& expected, double tolerance) {
      ASSERT_EQ(actual.size(), expected.size());
      for (size_t index = 0; index < expected.size(); ++index) {
         EXPECT_NEAR(actual[index], expected[index], tolerance) << "at index " << index;
      }
   }

   void expect_solution(const Result& result, const Solution& expected, double tolerance = 1e-6) {
      ASSERT_EQ(result.optimization_status, OptimizationStatus::SUCCESS);
      ASSERT_EQ(result.solution_status, SolutionStatus::FEASIBLE_KKT_POINT);
      EXPECT_NEAR(result.solution_objective, expected.objective, tolerance);
      expect_near(result.primal_solution, expected.primals, tolerance);
      expect_near(result.constraint_dual_solution, expected.constraint_multipliers, tolerance);
      expect_near(result.lower_bound_dual_solution, expected.lower_bound_multipliers, tolerance);
      expect_near(result.upper_bound_dual_solution, expected.upper_bound_multipliers, tolerance);
   }

   // stationarity ∇f - J^T y - z_L - z_U of the original problem (fixed variables included) at the returned solution
   std::vector<double> stationarity(const SparseTestProblem& problem, const Result& result) {
      Vector<double> gradient(problem.number_variables), jacobian(problem.jacobian_rows.size());
      problem.gradient(result.primal_solution, gradient);
      problem.jacobian(result.primal_solution, jacobian.data());
      std::vector<double> residual(problem.number_variables);
      for (size_t variable_index = 0; variable_index < problem.number_variables; ++variable_index) {
         residual[variable_index] = gradient[variable_index] - result.lower_bound_dual_solution[variable_index] -
            result.upper_bound_dual_solution[variable_index];
      }
      for (size_t nonzero_index = 0; nonzero_index < jacobian.size(); ++nonzero_index) {
         residual[static_cast<size_t>(problem.jacobian_columns[nonzero_index])] -= jacobian[nonzero_index] *
            result.constraint_dual_solution[static_cast<size_t>(problem.jacobian_rows[nonzero_index])];
      }
      return residual;
   }

   class Preprocessing: public ::testing::TestWithParam<std::string> {
   protected:
      void SetUp() override {
         if (!is_preset_available(GetParam())) {
            GTEST_SKIP() << "no subproblem solver available for the preset " << GetParam();
         }
      }
   };
} // namespace

TEST_P(Preprocessing, SolvesModelWithFixedVariables) {
   const SparseTestProblem problem = fixed_variable_problem();
   const SparseTestModel model(fixed_variable_problem());
   const Result result = solve(model, GetParam());
   expect_solution(result, {{1., 2., 2., 0.}, {-2., 1., 0., 0., 0.}, {0., 5., 0., 0.}, {0., 0., 0., 0.}, 8.});
   // the KKT conditions of the original problem hold at the returned solution
   for (double residual: stationarity(problem, result)) {
      EXPECT_NEAR(residual, 0., 1e-6);
   }
   EXPECT_LE(result.solution_primal_feasibility, 1e-8);
}

TEST_P(Preprocessing, SolvesModelWhoseVariablesAreAllFixed) {
   const SparseTestModel model(all_fixed_problem(-Inf));
   expect_solution(solve(model, GetParam()), {{1., 2.}, {0.}, {2., 1.}, {0., 0.}, 3.});
}

TEST_P(Preprocessing, ReportsInfeasibleModelWhoseVariablesAreAllFixed) {
   // x0 x1 = 2 >= 3 is violated: there is no feasible point
   const SparseTestModel model(all_fixed_problem(3.));
   const Result result = solve(model, GetParam());
   EXPECT_EQ(result.optimization_status, OptimizationStatus::SUCCESS);
   EXPECT_EQ(result.solution_status, SolutionStatus::INFEASIBLE_STATIONARY_POINT);
   expect_near(result.primal_solution, {1., 2.}, 0.);
}

// without iterations, no ingredient postprocesses the iterate: it must still be returned in the original space
TEST_P(Preprocessing, ReturnsInitialPointInOriginalSpaceWithoutIterations) {
   const SparseTestModel model(fixed_variable_problem());
   const Result result = solve(model, GetParam(), 0);
   EXPECT_EQ(result.optimization_status, OptimizationStatus::ITERATION_LIMIT);
   ASSERT_EQ(result.constraint_dual_solution.size(), size_t{5});
   ASSERT_EQ(result.constraint_values.size(), size_t{5});
   // initial point, with the fixed variable x1 at its value
   expect_near(result.primal_solution, {0., 2., 1., 0.5}, 0.);
}

INSTANTIATE_TEST_SUITE_P(Presets, Preprocessing, ::testing::Values("filtersqp", "ipopt"));
