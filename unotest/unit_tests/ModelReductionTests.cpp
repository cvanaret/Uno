// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <vector>
#include <gtest/gtest.h>
#include "../ReductionTestProblems.hpp"
#include "model/ModelReduction.hpp"
#include "model/ReducedModel.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Logger.hpp"

using namespace uno;

namespace {
   constexpr double feasibility_tolerance = 1e-9;

   template <typename Array>
   std::vector<double> to_vector(const Array& array, size_t size) {
      std::vector<double> vector(size);
      for (size_t index = 0; index < size; ++index) {
         vector[index] = array[index];
      }
      return vector;
   }

   void expect_near(const std::vector<double>& actual, const std::vector<double>& expected) {
      ASSERT_EQ(actual.size(), expected.size());
      for (size_t index = 0; index < actual.size(); ++index) {
         EXPECT_NEAR(actual[index], expected[index], 1e-12) << "at index " << index;
      }
   }
} // namespace

TEST(ModelReduction, EliminatesFixedVariablesAndConstraintsOfFixedVariablesOnly) {
   const SparseTestModel model(fixed_variable_problem());
   const ModelReduction reduction = ModelReduction::detect(model, feasibility_tolerance);
   EXPECT_FALSE(reduction.is_trivial());
   EXPECT_EQ(reduction.fixed_variables, (std::vector<size_t>{1}));
   EXPECT_EQ(reduction.free_variables, (std::vector<size_t>{0, 2, 3}));
   // c0 = x0 + x1 depends on a free variable: it is kept; c4 = x1^2 = 4 is satisfied (the initial point is x1 = 0)
   EXPECT_EQ(reduction.kept_constraints, (std::vector<size_t>{0, 1, 2, 3}));
   EXPECT_EQ(reduction.removed_constraints, (std::vector<size_t>{4}));
}

TEST(ModelReduction, KeepsViolatedConstraintsOfFixedVariablesOnly) {
   // x0 x1 = 2 >= 3 is violated: the solver reports the infeasibility
   const SparseTestModel model(all_fixed_problem(3.));
   const ModelReduction reduction = ModelReduction::detect(model, feasibility_tolerance);
   EXPECT_TRUE(reduction.free_variables.empty());
   EXPECT_EQ(reduction.kept_constraints, (std::vector<size_t>{0}));
   EXPECT_TRUE(reduction.removed_constraints.empty());
}

TEST(ModelReduction, KeepsConstraintsViolatedBeyondFeasibilityTolerance) {
   // x0 x1 = 2: a violation of 1e-12 is roundoff (tolerance 1e-12 max(1, |c|) = 2e-12), one of 2e-9 is not
   const SparseTestModel roundoff_model(all_fixed_problem(2. + 1e-12));
   EXPECT_EQ(ModelReduction::detect(roundoff_model, feasibility_tolerance).removed_constraints, (std::vector<size_t>{0}));
   const SparseTestModel violated_model(all_fixed_problem(2. + 2e-9));
   EXPECT_EQ(ModelReduction::detect(violated_model, feasibility_tolerance).kept_constraints, (std::vector<size_t>{0}));
   // 1e6 x0 x1 = 2e6 <= 2e6 - 1e-7: the roundoff tolerance (2e-6) is capped by the feasibility tolerance
   SparseTestProblem problem = all_fixed_problem(-Inf);
   problem.constraints = [](const Vector<double>& x, Vector<double>& c) { c[0] = 1e6*x[0]*x[1]; };
   problem.constraints_upper_bounds = {2e6 - 1e-7};
   const SparseTestModel large_model(std::move(problem));
   EXPECT_EQ(ModelReduction::detect(large_model, feasibility_tolerance).kept_constraints, (std::vector<size_t>{0}));
}

TEST(ModelReduction, IsTrivialWithoutFixedVariables) {
   SparseTestProblem problem = fixed_variable_problem();
   problem.variables_lower_bounds[1] = -Inf;
   const SparseTestModel model(std::move(problem));
   const ModelReduction reduction = ModelReduction::detect(model, feasibility_tolerance);
   EXPECT_TRUE(reduction.is_trivial());
   EXPECT_EQ(reduction.kept_constraints.size(), size_t{5});
}

TEST(ReducedModel, ReducesDimensionsSparsityAndEvaluations) {
   Logger::level = SILENT;
   const SparseTestModel model(fixed_variable_problem());
   const ReducedModel reduced_model(model, ModelReduction::detect(model, feasibility_tolerance));
   ASSERT_EQ(reduced_model.number_variables, size_t{3}); // x0, x2, x3
   ASSERT_EQ(reduced_model.number_constraints, size_t{4}); // c0, c1, c2, c3
   EXPECT_EQ(reduced_model.get_variables_lower_bounds(), (std::vector<double>{-Inf, 0., -Inf}));
   EXPECT_EQ(reduced_model.get_variables_upper_bounds(), (std::vector<double>{Inf, 10., Inf}));
   EXPECT_EQ(reduced_model.get_constraints_lower_bounds(), (std::vector<double>{-Inf, -4., 1., -Inf}));
   EXPECT_EQ(reduced_model.get_constraints_upper_bounds(), (std::vector<double>{3., Inf, 1., 4.}));
   EXPECT_EQ(reduced_model.get_equality_constraints(), (std::vector<size_t>{2}));
   EXPECT_EQ(reduced_model.get_inequality_constraints(), (std::vector<size_t>{0, 1, 3}));
   EXPECT_EQ(reduced_model.get_linear_constraints().size(), size_t{3});
   EXPECT_EQ(*reduced_model.get_nonlinear_constraints().begin(), size_t{3});

   // Jacobian: c0 w.r.t. x0, c1 w.r.t. x2, c2 w.r.t. x0 and x3, c3 w.r.t. x3
   ASSERT_EQ(reduced_model.number_jacobian_nonzeros(), size_t{5});
   expect_near(to_vector(reduced_model.get_jacobian_row_indices(), 5), {0., 1., 2., 2., 3.});
   expect_near(to_vector(reduced_model.get_jacobian_column_indices(), 5), {0., 1., 0., 2., 2.});
   // Hessian: the diagonal terms of x0, x2, x3 (1-based indexing)
   ASSERT_EQ(reduced_model.number_hessian_nonzeros(), size_t{3});
   std::vector<uno_int> rows(3), columns(3);
   reduced_model.compute_hessian_sparsity(view(rows.data(), 3), view(columns.data(), 3), 1);
   EXPECT_EQ(rows, (std::vector<uno_int>{1, 2, 3}));
   EXPECT_EQ(columns, (std::vector<uno_int>{1, 2, 3}));

   // evaluations at (x0, x2, x3) = (0.5, 1.5, -0.25), with x1 = 2
   Vector<double> x(3);
   x[0] = 0.5; x[1] = 1.5; x[2] = -0.25;
   EXPECT_DOUBLE_EQ(reduced_model.evaluate_objective(x), 6.25 + 1. + 2.25 + 0.0625 + 1.);
   Vector<double> constraints(4), gradient(3), jacobian(5), multipliers(4), hessian(3);
   reduced_model.evaluate_constraints(x, constraints);
   expect_near(to_vector(constraints, 4), {2.5, -3., 0.25, 0.0625});
   reduced_model.evaluate_objective_gradient(x, gradient);
   expect_near(to_vector(gradient, 3), {-3., -3., -0.5});
   reduced_model.evaluate_jacobian(x, jacobian.data());
   expect_near(to_vector(jacobian, 5), {1., -2., 1., 1., -0.5});
   multipliers[3] = 3.;
   reduced_model.evaluate_lagrangian_hessian(x, 1., multipliers, hessian.view());
   expect_near(to_vector(hessian, 3), {2., 2., -4.});
   Vector<double> initial_point(3);
   reduced_model.initial_primal_point(initial_point);
   expect_near(to_vector(initial_point, 3), {0., 1., 0.5});
}

TEST(ReducedModel, RecoversOriginalPrimalDualSolution) {
   Logger::level = SILENT;
   const SparseTestModel model(fixed_variable_problem());
   const ReducedModel reduced_model(model, ModelReduction::detect(model, feasibility_tolerance));
   // reduced KKT point: c0 (x0 + 2 <= 3) and c1 (x2 <= 2) are active
   Iterate iterate(reduced_model.number_variables, reduced_model.number_constraints);
   iterate.primals[0] = 1.; iterate.primals[1] = 2.; iterate.primals[2] = 0.;
   iterate.multipliers.constraints[0] = -2.; iterate.multipliers.constraints[1] = 1.;
   Evaluations evaluations(reduced_model);
   reduced_model.postprocess_solution(iterate, evaluations);
   expect_near(to_vector(iterate.primals, 4), {1., 2., 2., 0.});
   expect_near(to_vector(iterate.multipliers.constraints, 5), {-2., 1., 0., 0., 0.});
   // the bound multiplier of x1 absorbs the stationarity: ∇f_1 - J_01 y_0 = 3 + 2
   expect_near(to_vector(iterate.multipliers.lower_bounds, 4), {0., 5., 0., 0.});
   expect_near(to_vector(iterate.multipliers.upper_bounds, 4), {0., 0., 0., 0.});
}
