// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// The preprocessing (elimination of the fixed variables, bound constraints turned into variable bounds, detection of
// the linear constraints) is transparent: the KKT point is returned in the space of the original model

#include <algorithm>
#include <cmath>
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

   Result solve(const Model& model, const std::string& preset, bool detect_bound_constraints = true,
         bool detect_linear_constraints = true) {
      Options options;
      DefaultOptions::load(options);
      Presets::set(options, preset);
      options.set_bool("detect_bound_constraints", detect_bound_constraints);
      options.set_bool("detect_linear_constraints", detect_linear_constraints);
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

   class Preprocessing: public ::testing::TestWithParam<std::string> {
   protected:
      void SetUp() override {
         if (!is_preset_available(GetParam())) {
            GTEST_SKIP() << "no subproblem solver available for the preset " << GetParam();
         }
      }
   };
} // namespace

TEST_P(Preprocessing, SolvesModelWithFixedVariablesAndBoundConstraints) {
   const SparseTestModel model(mixed_problem());
   const Solution expected{{1., 2., 2., 0.}, {-2., 1., 0., 0., 0.}, {0., 5., 0., 0.}, {0., 0., 0., 0.}, 8.};
   expect_solution(solve(model, GetParam()), expected);
   // same KKT point without the optional detections
   expect_solution(solve(model, GetParam(), false, false), expected);
}

TEST_P(Preprocessing, SolvesModelWithVariablesFixedByBoundConstraints) {
   const SparseTestModel model(chained_problem());
   const Solution expected{{1., 2., 3.}, {6., -2.}, {0., 0., 0.}, {-4., 0., 0.}, 12.};
   expect_solution(solve(model, GetParam()), expected);
   expect_solution(solve(model, GetParam(), false, false), expected);
}

TEST_P(Preprocessing, SolvesModelWithUndeclaredLinearConstraints) {
   const SparseTestModel model(undeclared_linear_problem());
   const Solution expected{{2., -1., 2.}, {0., 0., -1.}, {0., 0., 2.}, {0., 0., 0.}, 1.};
   expect_solution(solve(model, GetParam()), expected);
   expect_solution(solve(model, GetParam(), true, false), expected);
}

INSTANTIATE_TEST_SUITE_P(Presets, Preprocessing, ::testing::Values("filtersqp", "ipopt"));

// the interior-point method relaxes the bound x <= 1 induced by the badly scaled constraint 1e10 x <= 1e10: the solution
// is projected onto the induced bound, otherwise the original constraint would be violated by up to 1e10 x relaxation
TEST(PreprocessingBoundRelaxation, BadlyScaledBoundConstraintIsSatisfied) {
   if (!is_preset_available("ipopt")) {
      GTEST_SKIP() << "no linear solver available";
   }
   constexpr double coefficient = 1e10;
   SparseTestProblem problem;
   problem.name = "badly_scaled_bound_constraint";
   problem.number_variables = 1;
   problem.number_constraints = 1;
   problem.objective = [](const Vector<double>& x) { return (x[0] - 2.)*(x[0] - 2.); };
   problem.gradient = [](const Vector<double>& x, Vector<double>& g) { g[0] = 2.*(x[0] - 2.); };
   problem.constraints = [](const Vector<double>& x, Vector<double>& c) { c[0] = coefficient * x[0]; };
   problem.jacobian_rows = {0};
   problem.jacobian_columns = {0};
   problem.jacobian = [](const Vector<double>&, double* J) { J[0] = coefficient; };
   problem.hessian_rows = {0};
   problem.hessian_columns = {0};
   problem.hessian = [](const Vector<double>&, double rho, const Vector<double>&, double* H) { H[0] = 2.*rho; };
   problem.variables_lower_bounds = {-Inf};
   problem.variables_upper_bounds = {Inf};
   problem.constraints_lower_bounds = {-Inf};
   problem.constraints_upper_bounds = {coefficient};
   problem.initial_point = {0.};
   problem.linear_constraints = {0};
   const SparseTestModel model(std::move(problem));

   const Result result = solve(model, "ipopt");
   ASSERT_EQ(result.optimization_status, OptimizationStatus::SUCCESS);
   EXPECT_LE(result.primal_solution[0], 1.);
   EXPECT_LE(result.constraint_values[0], coefficient);
   EXPECT_NEAR(result.constraint_dual_solution[0], -2. / coefficient, 1e-12);
}

// all the variables are eliminated: x0 = 1 (fixed), x0 + x1 = 3 (x1 = 2), x0 x1 <= 5 (kept, satisfied)
// min x0^2 + x1 => f = 3, y = (1, 0), z = (0, 0) (x1: 1 - y0 = 0; x0: 2 x0 - y0 - x1 y1 = 1)
TEST_P(Preprocessing, SolvesModelWhoseVariablesAreAllEliminated) {
   SparseTestProblem problem;
   problem.name = "all_eliminated";
   problem.number_variables = 2;
   problem.number_constraints = 2;
   problem.objective = [](const Vector<double>& x) { return x[0]*x[0] + x[1]; };
   problem.gradient = [](const Vector<double>& x, Vector<double>& g) { g[0] = 2.*x[0]; g[1] = 1.; };
   problem.constraints = [](const Vector<double>& x, Vector<double>& c) { c[0] = x[0] + x[1]; c[1] = x[0]*x[1]; };
   problem.jacobian_rows = {0, 0, 1, 1};
   problem.jacobian_columns = {0, 1, 0, 1};
   problem.jacobian = [](const Vector<double>& x, double* J) { J[0] = 1.; J[1] = 1.; J[2] = x[1]; J[3] = x[0]; };
   problem.hessian_rows = {0, 1};
   problem.hessian_columns = {0, 0};
   problem.hessian = [](const Vector<double>&, double rho, const Vector<double>& y, double* H) { H[0] = 2.*rho; H[1] = -y[1]; };
   problem.variables_lower_bounds = {1., -Inf};
   problem.variables_upper_bounds = {1., Inf};
   problem.constraints_lower_bounds = {3., -Inf};
   problem.constraints_upper_bounds = {3., 5.};
   problem.initial_point = {0., 0.};
   problem.linear_constraints = {0};
   const SparseTestModel model(std::move(problem));
   const Result result = solve(model, GetParam());
   expect_solution(result, {{1., 2.}, {1., 0.}, {1., 0.}, {0., 0.}, 3.});
}

// without iterations, no ingredient postprocesses the iterate: it must still be returned in the original space
TEST(PreprocessingWithoutIterations, InitialPointIsReturnedInOriginalSpace) {
   const SparseTestModel model(mixed_problem());
   Options options;
   DefaultOptions::load(options);
   Presets::set(options, "ipopt");
   options.set_integer("max_iterations", 0);
   Logger::level = SILENT;
   Uno uno;
   const Result result = uno.solve(model, options);
   EXPECT_EQ(result.optimization_status, OptimizationStatus::ITERATION_LIMIT);
   ASSERT_EQ(result.primal_solution.size(), size_t{4});
   ASSERT_EQ(result.constraint_dual_solution.size(), size_t{5});
   ASSERT_EQ(result.lower_bound_dual_solution.size(), size_t{4});
   ASSERT_EQ(result.constraint_values.size(), size_t{5});
   // initial point, with the fixed variable x1 at its value
   expect_near(result.primal_solution, {0., 2., 1., 0.5}, 0.);
}
