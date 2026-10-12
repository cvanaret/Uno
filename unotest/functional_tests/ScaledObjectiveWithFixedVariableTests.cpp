// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Regression test: the objective and constraint values reported in the Result must be those of the original
// (unscaled) model at the returned primal solution.
// Bug: with function scaling + interior-point method + a fixed variable, PrimalDualInteriorPointProblem::postprocess_solution
// snaps the fixed variable and resets the evaluations; ScaledModel::postprocess_solution then skips the unscaling of
// the (already evaluated) objective, so Result::solution_objective is the *scaled* objective sigma_f * f(x*).

#include <cmath>
#include <vector>
#include <gtest/gtest.h>
#include "Uno.hpp"
#include "model/Model.hpp"
#include "optimization/Result.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"
#include "symbolic/IntegerRange.hpp"
#include "tools/Logger.hpp"

using namespace uno;

namespace {
   // min  1000 [(x0 - 1)^2 + (x1 - 2)^2 + x0 x2 + x1 x2]
   // s.t. x0^2 + x1^2 + x2 <= 4
   //      x2 = 0 (fixed variable)
   // Large objective gradient (> function_scaling_threshold = 100) => the objective is scaled.
   // The coupling with x2 makes the Newton step on the fixed variable subject to roundoff, so that it gets snapped.
   class FixedVariableModel: public Model {
   public:
      FixedVariableModel(): Model("fixed_variable", 3, 1, 1., -1., 0) { }

      ProblemType get_problem_type() const override { return ProblemType::NONLINEAR; }
      bool has_jacobian_operator() const override { return false; }
      bool has_jacobian_transposed_operator() const override { return false; }
      bool has_hessian_operator() const override { return false; }
      bool has_hessian_matrix() const override { return true; }

      double evaluate_objective(const Vector<double>& x) const override {
         ++this->number_objective_evaluations;
         return 1000. * (std::pow(x[0] - 1., 2) + std::pow(x[1] - 2., 2) + x[0] * x[2] + x[1] * x[2]);
      }

      void evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const override {
         ++this->number_constraint_evaluations;
         constraints[0] = x[0] * x[0] + x[1] * x[1] + x[2];
      }

      void evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const override {
         ++this->number_gradient_evaluations;
         gradient[0] = 1000. * (2. * (x[0] - 1.) + x[2]);
         gradient[1] = 1000. * (2. * (x[1] - 2.) + x[2]);
         gradient[2] = 1000. * (x[0] + x[1]);
      }

      View<const uno_int> get_jacobian_row_indices() const override {
         return view(this->jacobian_rows.data(), this->jacobian_rows.size());
      }

      View<const uno_int> get_jacobian_column_indices() const override {
         return view(this->jacobian_columns.data(), this->jacobian_columns.size());
      }

      void evaluate_jacobian(const Vector<double>& x, double* jacobian_values) const override {
         ++this->number_jacobian_evaluations;
         jacobian_values[0] = 2. * x[0];
         jacobian_values[1] = 2. * x[1];
         jacobian_values[2] = 1.;
      }

      // lower triangle: (0,0), (1,1), (2,0), (2,1)
      void compute_hessian_sparsity(View<uno_int> row_indices, View<uno_int> column_indices, uno_int solver_indexing) const override {
         const uno_int rows[] = {0, 1, 2, 2};
         const uno_int columns[] = {0, 1, 0, 1};
         for (size_t index = 0; index < 4; ++index) {
            row_indices[index] = rows[index] + solver_indexing;
            column_indices[index] = columns[index] + solver_indexing;
         }
      }

      // Lagrangian rho f(x) - y^T c(x)
      void evaluate_lagrangian_hessian(const Vector<double>& /*x*/, double objective_multiplier, const Vector<double>& multipliers,
            View<double> hessian_values) const override {
         ++this->number_hessian_evaluations;
         hessian_values[0] = objective_multiplier * 2000. - multipliers[0] * 2.;
         hessian_values[1] = objective_multiplier * 2000. - multipliers[0] * 2.;
         hessian_values[2] = objective_multiplier * 1000.;
         hessian_values[3] = objective_multiplier * 1000.;
      }

      void compute_jacobian_vector_product(const double*, const double*, double*) const override {
         throw std::runtime_error("not implemented");
      }
      void compute_jacobian_transposed_vector_product(const double*, const double*, double*) const override {
         throw std::runtime_error("not implemented");
      }
      void compute_hessian_vector_product(View<const double>, View<const double>, double, const Vector<double>&,
            View<double>) const override {
         throw std::runtime_error("not implemented");
      }

      const std::vector<double>& get_variables_lower_bounds() const override { return this->variables_lower_bounds; }
      const std::vector<double>& get_variables_upper_bounds() const override { return this->variables_upper_bounds; }
      const std::vector<double>& get_constraints_lower_bounds() const override { return this->constraints_lower_bounds; }
      const std::vector<double>& get_constraints_upper_bounds() const override { return this->constraints_upper_bounds; }
      const std::vector<size_t>& get_equality_constraints() const override { return this->equality_constraints; }
      const std::vector<size_t>& get_inequality_constraints() const override { return this->inequality_constraints; }
      const Collection<size_t>& get_linear_constraints() const override { return this->linear_constraints; }
      const Collection<size_t>& get_nonlinear_constraints() const override { return this->nonlinear_constraints; }

      void initial_primal_point(Vector<double>& x) const override {
         x[0] = 0.5;
         x[1] = 0.5;
         x[2] = 0.3; // not on its fixed value
      }
      void initial_dual_point(Vector<double>& multipliers) const override { multipliers.fill(0.); }
      void postprocess_solution(Iterate& /*iterate*/) const override { }

      size_t number_jacobian_nonzeros() const override { return 3; }
      size_t number_hessian_nonzeros() const override { return 4; }

      size_t number_model_objective_evaluations() const override { return this->number_objective_evaluations; }
      size_t number_model_constraints_evaluations() const override { return this->number_constraint_evaluations; }
      size_t number_model_objective_gradient_evaluations() const override { return this->number_gradient_evaluations; }
      size_t number_model_jacobian_evaluations() const override { return this->number_jacobian_evaluations; }
      size_t number_model_hessian_evaluations() const override { return this->number_hessian_evaluations; }
      void reset_number_evaluations() const override {
         this->number_objective_evaluations = this->number_constraint_evaluations = this->number_gradient_evaluations =
            this->number_jacobian_evaluations = this->number_hessian_evaluations = 0;
      }

   private:
      const std::vector<uno_int> jacobian_rows{0, 0, 0};
      const std::vector<uno_int> jacobian_columns{0, 1, 2};
      const std::vector<double> variables_lower_bounds{-INFINITY, -INFINITY, 0.};
      const std::vector<double> variables_upper_bounds{INFINITY, INFINITY, 0.};
      const std::vector<double> constraints_lower_bounds{-INFINITY};
      const std::vector<double> constraints_upper_bounds{4.};
      const std::vector<size_t> equality_constraints{};
      const std::vector<size_t> inequality_constraints{0};
      const IntegerRange linear_constraints{0};
      const IntegerRange nonlinear_constraints{1};
      mutable size_t number_objective_evaluations{0};
      mutable size_t number_constraint_evaluations{0};
      mutable size_t number_gradient_evaluations{0};
      mutable size_t number_jacobian_evaluations{0};
      mutable size_t number_hessian_evaluations{0};
   };
} // namespace

TEST(ScaledObjectiveWithFixedVariable, ReportedObjectiveAndConstraintsAreUnscaled) {
   const FixedVariableModel model;
   Options options;
   DefaultOptions::load(options);
   Presets::set(options, "ipopt"); // interior point + log barrier + use_function_scaling = true
   Logger::level = SILENT;

   Uno uno;
   const Result result = uno.solve(model, options);
   ASSERT_EQ(result.optimization_status, OptimizationStatus::SUCCESS);
   ASSERT_EQ(result.solution_status, SolutionStatus::FEASIBLE_KKT_POINT);

   // the fixed variable is returned on its value
   EXPECT_EQ(result.primal_solution[2], 0.);

   // reported values must be those of the original model at the reported solution
   const double objective = model.evaluate_objective(result.primal_solution);
   EXPECT_NEAR(result.solution_objective, objective, 1e-8 * std::max(1., std::abs(objective)));

   Vector<double> constraints(1);
   model.evaluate_constraints(result.primal_solution, constraints);
   EXPECT_NEAR(result.constraint_values[0], constraints[0], 1e-10);

   // sanity check against the analytical solution: x = 2 (1, 2)/sqrt(5), f = 1000 (sqrt(5) - 2)^2
   const double expected_objective = 1000. * std::pow(std::sqrt(5.) - 2., 2);
   EXPECT_NEAR(result.solution_objective, expected_objective, 1e-6 * expected_objective);
}
