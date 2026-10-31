// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <vector>
#include <gtest/gtest.h>
#include "../ReductionTestProblems.hpp"
#include "model/ModelReduction.hpp"
#include "model/ReducedModel.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"

using namespace uno;

namespace {
   Options default_options(bool detect_bound_constraints = true, bool detect_linear_constraints = true) {
      Options options;
      DefaultOptions::load(options);
      options.set_bool("detect_bound_constraints", detect_bound_constraints);
      options.set_bool("detect_linear_constraints", detect_linear_constraints);
      Logger::level = SILENT;
      return options;
   }

   template <typename Array>
   std::vector<double> to_vector(const Array& array, size_t size) {
      std::vector<double> vector(size);
      for (size_t index = 0; index < size; ++index) {
         vector[index] = array[index];
      }
      return vector;
   }

   // the solution of the reduced problem, given by its primals and its multipliers
   Iterate reduced_iterate(const Model& reduced_model, const std::vector<double>& primals, const std::vector<double>& constraint_multipliers,
         const std::vector<double>& upper_bound_multipliers) {
      Iterate iterate(reduced_model.number_variables, reduced_model.number_constraints);
      for (size_t index = 0; index < primals.size(); ++index) {
         iterate.primals[index] = primals[index];
         iterate.multipliers.upper_bounds[index] = upper_bound_multipliers[index];
      }
      for (size_t index = 0; index < constraint_multipliers.size(); ++index) {
         iterate.multipliers.constraints[index] = constraint_multipliers[index];
      }
      return iterate;
   }

   void expect_near(const std::vector<double>& actual, const std::vector<double>& expected, double tolerance = 1e-12) {
      ASSERT_EQ(actual.size(), expected.size());
      for (size_t index = 0; index < actual.size(); ++index) {
         EXPECT_NEAR(actual[index], expected[index], tolerance) << "at index " << index;
      }
   }
} // namespace

TEST(ModelReduction, DetectsFixedVariablesAndBoundConstraints) {
   const SparseTestModel model(mixed_problem());
   const ModelReduction reduction = ModelReduction::detect(model, default_options());
   EXPECT_EQ(reduction.fixed_variables, (std::vector<size_t>{1}));
   EXPECT_EQ(reduction.elimination_rounds, (std::vector<size_t>{0}));
   EXPECT_EQ(reduction.free_variables, (std::vector<size_t>{0, 2, 3}));
   EXPECT_EQ(reduction.kept_constraints, (std::vector<size_t>{2, 3}));
   EXPECT_TRUE(reduction.detected_linear_constraints.empty()); // x3^2 is nonlinear
   ASSERT_EQ(reduction.bound_constraints.size(), size_t{3});
   const std::vector<std::pair<size_t, size_t>> expected{{0, 0}, {1, 2}, {4, 2}};
   const std::vector<double> coefficients{1., -2., 3.};
   for (size_t index = 0; index < 3; ++index) {
      EXPECT_EQ(reduction.bound_constraints[index].constraint_index, expected[index].first);
      EXPECT_EQ(reduction.bound_constraints[index].variable_index, expected[index].second);
      EXPECT_EQ(reduction.bound_constraints[index].coefficient, coefficients[index]);
   }
   // x0 <= 1 from c0 (with x1 at its value, although the initial point is x1 = 0)
   EXPECT_EQ(reduction.variables_upper_bounds[0], 1.);
   EXPECT_EQ(reduction.upper_bound_origins[0], std::optional<size_t>{0});
   // x2 <= 2 from c1; x2 >= -10 from c4 is looser than x2 >= 0
   EXPECT_EQ(reduction.variables_lower_bounds[2], 0.);
   EXPECT_FALSE(reduction.lower_bound_origins[2].has_value());
   EXPECT_EQ(reduction.variables_upper_bounds[2], 2.);
   EXPECT_EQ(reduction.upper_bound_origins[2], std::optional<size_t>{1});
}

TEST(ModelReduction, EliminatesVariablesFixedByBoundConstraints) {
   const SparseTestModel model(chained_problem());
   const ModelReduction reduction = ModelReduction::detect(model, default_options());
   EXPECT_EQ(reduction.fixed_variables, (std::vector<size_t>{0, 1}));
   EXPECT_EQ(reduction.elimination_rounds, (std::vector<size_t>{0, 1}));
   EXPECT_EQ(reduction.free_variables, (std::vector<size_t>{2}));
   EXPECT_TRUE(reduction.kept_constraints.empty());
   EXPECT_EQ(reduction.variables_lower_bounds[1], 2.);
   EXPECT_EQ(reduction.variables_upper_bounds[1], 2.);
   EXPECT_EQ(reduction.variables_upper_bounds[2], 3.);
   EXPECT_EQ(reduction.upper_bound_origins[2], std::optional<size_t>{1});
}

TEST(ModelReduction, DetectsUndeclaredLinearConstraints) {
   const SparseTestModel model(undeclared_linear_problem());
   const ModelReduction reduction = ModelReduction::detect(model, default_options());
   // x0 x2 is linear once x2 is fixed
   EXPECT_EQ(reduction.detected_linear_constraints, (std::vector<size_t>{0, 2}));
   ASSERT_EQ(reduction.bound_constraints.size(), size_t{1});
   EXPECT_EQ(reduction.bound_constraints[0].constraint_index, size_t{2});
   EXPECT_EQ(reduction.bound_constraints[0].coefficient, 2.);
   EXPECT_EQ(reduction.variables_upper_bounds[0], 2.);
   EXPECT_EQ(reduction.kept_constraints, (std::vector<size_t>{0, 1}));
}

TEST(ModelReduction, DetectionsCanBeDisabled) {
   const SparseTestModel mixed_model(mixed_problem());
   const ModelReduction reduction = ModelReduction::detect(mixed_model, default_options(false, false));
   // the fixed variables are always eliminated
   EXPECT_EQ(reduction.fixed_variables, (std::vector<size_t>{1}));
   EXPECT_TRUE(reduction.bound_constraints.empty());
   EXPECT_EQ(reduction.kept_constraints.size(), size_t{5});

   const SparseTestModel undeclared_linear_model(undeclared_linear_problem());
   const ModelReduction undeclared_reduction = ModelReduction::detect(undeclared_linear_model, default_options(true, false));
   EXPECT_TRUE(undeclared_reduction.detected_linear_constraints.empty());
   EXPECT_TRUE(undeclared_reduction.bound_constraints.empty());
}

TEST(ModelReduction, IsTrivialWithoutFixedVariablesAndBoundConstraints) {
   SparseTestProblem problem = mixed_problem();
   problem.variables_lower_bounds[1] = -Inf;
   const SparseTestModel model(std::move(problem));
   EXPECT_TRUE(ModelReduction::detect(model, default_options(false, true)).is_trivial());
   EXPECT_FALSE(ModelReduction::detect(model, default_options(true, true)).is_trivial());
}

TEST(ReducedModel, ReducesDimensionsSparsityAndEvaluations) {
   const SparseTestModel model(mixed_problem());
   const ReducedModel reduced_model(model, ModelReduction::detect(model, default_options()));
   ASSERT_EQ(reduced_model.number_variables, size_t{3}); // x0, x2, x3
   ASSERT_EQ(reduced_model.number_constraints, size_t{2}); // c2, c3
   EXPECT_EQ(reduced_model.get_variables_lower_bounds(), (std::vector<double>{-Inf, 0., -Inf}));
   EXPECT_EQ(reduced_model.get_variables_upper_bounds(), (std::vector<double>{1., 2., Inf}));
   EXPECT_EQ(reduced_model.get_constraints_lower_bounds(), (std::vector<double>{1., -Inf}));
   EXPECT_EQ(reduced_model.get_constraints_upper_bounds(), (std::vector<double>{1., 4.}));
   EXPECT_EQ(reduced_model.get_equality_constraints(), (std::vector<size_t>{0}));
   EXPECT_EQ(reduced_model.get_inequality_constraints(), (std::vector<size_t>{1}));
   EXPECT_EQ(reduced_model.get_linear_constraints().size(), size_t{1});
   EXPECT_EQ(*reduced_model.get_nonlinear_constraints().begin(), size_t{1});

   // Jacobian: c2 w.r.t. x0 and x3, c3 w.r.t. x3
   ASSERT_EQ(reduced_model.number_jacobian_nonzeros(), size_t{3});
   expect_near(to_vector(reduced_model.get_jacobian_row_indices(), 3), {0., 0., 1.});
   expect_near(to_vector(reduced_model.get_jacobian_column_indices(), 3), {0., 2., 2.});
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
   Vector<double> constraints(2), gradient(3), jacobian(3), multipliers(2), hessian(3);
   reduced_model.evaluate_constraints(x, constraints);
   expect_near(to_vector(constraints, 2), {0.25, 0.0625});
   reduced_model.evaluate_objective_gradient(x, gradient);
   expect_near(to_vector(gradient, 3), {-3., -3., -0.5});
   reduced_model.evaluate_jacobian(x, jacobian.data());
   expect_near(to_vector(jacobian, 3), {1., 1., -0.5});
   multipliers[1] = 3.;
   reduced_model.evaluate_lagrangian_hessian(x, 1., multipliers, hessian.view());
   expect_near(to_vector(hessian, 3), {2., 2., -4.});
   Vector<double> initial_point(3);
   reduced_model.initial_primal_point(initial_point);
   expect_near(to_vector(initial_point, 3), {0., 1., 0.5});
}

TEST(ReducedModel, RecoversOriginalPrimalDualSolution) {
   const SparseTestModel model(mixed_problem());
   const ReducedModel reduced_model(model, ModelReduction::detect(model, default_options()));
   // reduced KKT point: the upper bounds x0 <= 1 and x2 <= 2 are active
   Iterate iterate = reduced_iterate(reduced_model, {1., 2., 0.}, {0., 0.}, {-2., -2., 0.});
   Evaluations evaluations(reduced_model);
   reduced_model.postprocess_solution(iterate, evaluations);
   expect_near(to_vector(iterate.primals, 4), {1., 2., 2., 0.});
   expect_near(to_vector(iterate.multipliers.constraints, 5), {-2., 1., 0., 0., 0.});
   expect_near(to_vector(iterate.multipliers.lower_bounds, 4), {0., 5., 0., 0.});
   expect_near(to_vector(iterate.multipliers.upper_bounds, 4), {0., 0., 0., 0.});
}

TEST(ReducedModel, RecoversMultipliersAcrossEliminationRounds) {
   const SparseTestModel model(chained_problem());
   const ReducedModel reduced_model(model, ModelReduction::detect(model, default_options()));
   ASSERT_EQ(reduced_model.number_variables, size_t{1});
   ASSERT_EQ(reduced_model.number_constraints, size_t{0});
   Iterate iterate = reduced_iterate(reduced_model, {3.}, {}, {-2.});
   Evaluations evaluations(reduced_model);
   reduced_model.postprocess_solution(iterate, evaluations);
   expect_near(to_vector(iterate.primals, 3), {1., 2., 3.});
   expect_near(to_vector(iterate.multipliers.constraints, 2), {6., -2.});
   expect_near(to_vector(iterate.multipliers.lower_bounds, 3), {0., 0., 0.});
   expect_near(to_vector(iterate.multipliers.upper_bounds, 3), {-4., 0., 0.});
}

TEST(ReducedModel, ProjectsOntoBoundsInducedByConstraints) {
   const SparseTestModel model(mixed_problem());
   const ReducedModel reduced_model(model, ModelReduction::detect(model, default_options()));
   // slight violations of the relaxed bounds x0 <= 1 and x2 <= 2
   Iterate iterate = reduced_iterate(reduced_model, {1. + 1e-7, 2. + 1e-7, 0.}, {0., 0.}, {0., 0., 0.});
   Evaluations evaluations(reduced_model);
   evaluations.is_objective_computed = true;
   reduced_model.postprocess_solution(iterate, evaluations);
   EXPECT_EQ(iterate.primals[0], 1.);
   EXPECT_EQ(iterate.primals[2], 2.);
   EXPECT_FALSE(evaluations.is_objective_computed);
}

// x0 = 1 and x1 = 2 fixed: c0 = x0 + x1 <= 4 is a satisfied constant (removed), c1 = x0 - x1 >= 0 is a violated constant
// (kept, the solver reports the infeasibility), c2 = x0 + x2 >= 0 is the bound x2 >= -1
TEST(ModelReduction, RemovesSatisfiedConstantConstraints) {
   SparseTestProblem problem;
   problem.name = "constant_constraints";
   problem.number_variables = 3;
   problem.number_constraints = 3;
   problem.objective = [](const Vector<double>& x) { return x[2]*x[2]; };
   problem.gradient = [](const Vector<double>& x, Vector<double>& g) { g[0] = 0.; g[1] = 0.; g[2] = 2.*x[2]; };
   problem.constraints = [](const Vector<double>& x, Vector<double>& c) { c[0] = x[0] + x[1]; c[1] = x[0] - x[1]; c[2] = x[0] + x[2]; };
   problem.jacobian_rows = {0, 0, 1, 1, 2, 2};
   problem.jacobian_columns = {0, 1, 0, 1, 0, 2};
   problem.jacobian = [](const Vector<double>&, double* J) { J[0] = 1.; J[1] = 1.; J[2] = 1.; J[3] = -1.; J[4] = 1.; J[5] = 1.; };
   problem.hessian_rows = {2};
   problem.hessian_columns = {2};
   problem.hessian = [](const Vector<double>&, double rho, const Vector<double>&, double* H) { H[0] = 2.*rho; };
   problem.variables_lower_bounds = {1., 2., -Inf};
   problem.variables_upper_bounds = {1., 2., Inf};
   problem.constraints_lower_bounds = {-Inf, 0., 0.};
   problem.constraints_upper_bounds = {4., Inf, Inf};
   problem.initial_point = {0., 0., 0.};
   problem.linear_constraints = {0, 1, 2};
   const SparseTestModel model(std::move(problem));
   const ModelReduction reduction = ModelReduction::detect(model, default_options());
   EXPECT_EQ(reduction.redundant_constraints, (std::vector<size_t>{0}));
   EXPECT_EQ(reduction.kept_constraints, (std::vector<size_t>{1}));
   ASSERT_EQ(reduction.bound_constraints.size(), size_t{1});
   EXPECT_EQ(reduction.variables_lower_bounds[2], -1.);
}

// 3 x0 >= 0.3 induces x0 >= 0.3/3, which differs from the bound x0 <= 0.1 by roundoff only: x0 is fixed at 0.1
TEST(ModelReduction, FixesVariableWithAlmostEqualInducedBounds) {
   const SparseTestModel model(almost_fixed_problem());
   const ModelReduction reduction = ModelReduction::detect(model, default_options());
   EXPECT_EQ(reduction.fixed_variables, (std::vector<size_t>{0}));
   EXPECT_EQ(reduction.variables_lower_bounds[0], 0.1);
   EXPECT_EQ(reduction.variables_upper_bounds[0], 0.1);
   EXPECT_TRUE(reduction.kept_constraints.empty());
}

// a satisfied nonlinear constraint of the fixed variables only would have an empty Jacobian row: it is removed
TEST(ModelReduction, RemovesSatisfiedConstraintsOfFixedVariablesOnly) {
   const SparseTestModel model(fixed_only_constraint_problem());
   const ModelReduction reduction = ModelReduction::detect(model, default_options(false, false));
   EXPECT_EQ(reduction.redundant_constraints, (std::vector<size_t>{0}));
   EXPECT_EQ(reduction.kept_constraints, (std::vector<size_t>{1}));
}
