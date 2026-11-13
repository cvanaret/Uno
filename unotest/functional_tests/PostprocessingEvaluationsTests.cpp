// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Regression test: the postprocessing must not call the user functions again when the final iterate is unchanged.
// An extra objective/constraint call after the last accepted iterate becomes the last record of drivers that log
// every evaluation (e.g. the pyOptSparse history), without the major-iteration data attached to it.

#include <cmath>
#include <vector>
#include <gtest/gtest.h>
#include "Uno.hpp"
#include "model/Model.hpp"
#include "optimization/Multipliers.hpp"
#include "optimization/Result.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"
#include "symbolic/IntegerRange.hpp"
#include "tools/Logger.hpp"
#include "tools/UserCallbacks.hpp"

using namespace uno;

namespace {
   class HS15Model: public Model {
   public:
      HS15Model(): Model("hs015", 2, 2, 1., -1., 0) { }

      ProblemType get_problem_type() const override { return ProblemType::NONLINEAR; }
      bool has_jacobian_operator() const override { return false; }
      bool has_jacobian_transposed_operator() const override { return false; }
      bool has_hessian_operator() const override { return false; }
      bool has_hessian_matrix() const override { return true; }

      double evaluate_objective(const Vector<double>& x) const override {
         ++this->number_objective_evaluations;
         return 100. * std::pow(x[1] - x[0] * x[0], 2) + std::pow(1. - x[0], 2);
      }

      void evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const override {
         ++this->number_constraint_evaluations;
         constraints[0] = x[0] * x[1];
         constraints[1] = x[0] + x[1] * x[1];
      }

      void evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const override {
         ++this->number_gradient_evaluations;
         gradient[0] = -400. * x[0] * (x[1] - x[0] * x[0]) - 2. * (1. - x[0]);
         gradient[1] = 200. * (x[1] - x[0] * x[0]);
      }

      View<const uno_int> get_jacobian_row_indices() const override {
         return view(this->jacobian_rows.data(), this->jacobian_rows.size());
      }

      View<const uno_int> get_jacobian_column_indices() const override {
         return view(this->jacobian_columns.data(), this->jacobian_columns.size());
      }

      void evaluate_jacobian(const Vector<double>& x, double* jacobian_values) const override {
         ++this->number_jacobian_evaluations;
         jacobian_values[0] = x[1];
         jacobian_values[1] = x[0];
         jacobian_values[2] = 1.;
         jacobian_values[3] = 2. * x[1];
      }

      // lower triangle: (0,0), (1,0), (1,1)
      void compute_hessian_sparsity(View<uno_int> row_indices, View<uno_int> column_indices, uno_int solver_indexing) const override {
         const uno_int rows[] = {0, 1, 1};
         const uno_int columns[] = {0, 0, 1};
         for (size_t index = 0; index < 3; ++index) {
            row_indices[index] = rows[index] + solver_indexing;
            column_indices[index] = columns[index] + solver_indexing;
         }
      }

      // Lagrangian rho f(x) - y^T c(x)
      void evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
            View<double> hessian_values) const override {
         ++this->number_hessian_evaluations;
         hessian_values[0] = objective_multiplier * (1200. * x[0] * x[0] - 400. * x[1] + 2.);
         hessian_values[1] = objective_multiplier * (-400. * x[0]) - multipliers[0];
         hessian_values[2] = objective_multiplier * 200. - multipliers[1] * 2.;
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
         x[0] = -2.;
         x[1] = 1.;
      }
      void initial_dual_point(Vector<double>& multipliers) const override { multipliers.fill(0.); }
      void postprocess_solution(Iterate& /*iterate*/, Evaluations& /*evaluations*/) const override { }

      size_t number_jacobian_nonzeros() const override { return 4; }
      size_t number_hessian_nonzeros() const override { return 3; }

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
      const std::vector<uno_int> jacobian_rows{0, 0, 1, 1};
      const std::vector<uno_int> jacobian_columns{0, 1, 0, 1};
      const std::vector<double> variables_lower_bounds{-INFINITY, -INFINITY};
      const std::vector<double> variables_upper_bounds{0.5, INFINITY};
      const std::vector<double> constraints_lower_bounds{1., 0.};
      const std::vector<double> constraints_upper_bounds{INFINITY, INFINITY};
      const std::vector<size_t> equality_constraints{};
      const std::vector<size_t> inequality_constraints{0, 1};
      const IntegerRange linear_constraints{0};
      const IntegerRange nonlinear_constraints{2};
   public:
      mutable size_t number_objective_evaluations{0};
      mutable size_t number_constraint_evaluations{0};
      mutable size_t number_gradient_evaluations{0};
      mutable size_t number_jacobian_evaluations{0};
      mutable size_t number_hessian_evaluations{0};
   };

   // records the number of user function evaluations when each acceptable iterate is notified
   class EvaluationCountCallbacks: public UserCallbacks {
   public:
      explicit EvaluationCountCallbacks(const HS15Model& model): model(model) { }

      void notify_acceptable_iterate(const Vector<double>& /*primals*/, const Multipliers& /*multipliers*/, double /*objective_multiplier*/,
            double /*primal_feasibility_residual*/, double /*stationarity_residual*/, double /*complementarity_residual*/) override {
         ++this->number_notifications;
         this->objective_evaluations_at_last_notification = this->model.number_objective_evaluations;
         this->constraint_evaluations_at_last_notification = this->model.number_constraint_evaluations;
      }
      bool termination(const Vector<double>& /*primals*/, const Multipliers& /*multipliers*/, double /*objective_multiplier*/,
            double /*primal_feasibility_residual*/, double /*stationarity_residual*/, double /*complementarity_residual*/) override {
         return false;
      }

      const HS15Model& model;
      size_t number_notifications{0};
      size_t objective_evaluations_at_last_notification{0};
      size_t constraint_evaluations_at_last_notification{0};
   };
} // namespace

TEST(Postprocessing, NoUserEvaluationAfterLastAcceptedIterate) {
   const HS15Model model;
   Options options;
   DefaultOptions::load(options);
   Presets::set(options, "ipopt"); // function scaling on
   Logger::level = SILENT;

   EvaluationCountCallbacks callbacks(model);
   Uno uno;
   const Result result = uno.solve(model, options, callbacks);
   ASSERT_EQ(result.optimization_status, OptimizationStatus::SUCCESS);
   ASSERT_EQ(result.solution_status, SolutionStatus::FEASIBLE_KKT_POINT);
   ASSERT_LT(0u, callbacks.number_notifications);

   EXPECT_EQ(model.number_objective_evaluations, callbacks.objective_evaluations_at_last_notification);
   EXPECT_EQ(model.number_constraint_evaluations, callbacks.constraint_evaluations_at_last_notification);

   // the cached values must still be reported in the units of the original model
   EXPECT_NEAR(result.solution_objective, model.evaluate_objective(result.primal_solution), 1e-10);
   Vector<double> constraints(2);
   model.evaluate_constraints(result.primal_solution, constraints);
   EXPECT_NEAR(result.constraint_values[0], constraints[0], 1e-10);
   EXPECT_NEAR(result.constraint_values[1], constraints[1], 1e-10);
   EXPECT_NEAR(result.solution_objective, 306.5, 1e-3); // HS15 optimum (0.5, 2)
}
