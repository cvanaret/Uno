// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_TEST_BOUNDACTIVEMODEL_H
#define UNO_TEST_BOUNDACTIVEMODEL_H

#include <cmath>
#include <stdexcept>
#include <vector>
#include "model/Model.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/IntegerRange.hpp"

namespace uno::test {
   // min (x0 + 1)^2 + (x1 - 2)^2  s.t.  x0 + x1 >= 1,  x0 >= 0,  -1 <= x1 <= 1
   // solution: x* = (0, 1) with both bounds active and the constraint active
   class BoundActiveModel: public Model {
   public:
      BoundActiveModel(): Model("bound_active", 2, 1, 1., -1., 0) { }

      ProblemType get_problem_type() const override { return ProblemType::NONLINEAR; }
      bool has_jacobian_operator() const override { return false; }
      bool has_jacobian_transposed_operator() const override { return false; }
      bool has_hessian_operator() const override { return false; }
      bool has_hessian_matrix() const override { return true; }

      double evaluate_objective(const Vector<double>& x) const override {
         ++this->number_objective_evaluations;
         return std::pow(x[0] + 1., 2) + std::pow(x[1] - 2., 2);
      }
      void evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const override {
         ++this->number_constraint_evaluations;
         constraints[0] = x[0] + x[1];
      }
      void evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const override {
         gradient[0] = 2. * (x[0] + 1.);
         gradient[1] = 2. * (x[1] - 2.);
      }
      View<const uno_int> get_jacobian_row_indices() const override { return view(this->jacobian_rows.data(), 2); }
      View<const uno_int> get_jacobian_column_indices() const override { return view(this->jacobian_columns.data(), 2); }
      void evaluate_jacobian(const Vector<double>& /*x*/, double* jacobian_values) const override {
         jacobian_values[0] = 1.;
         jacobian_values[1] = 1.;
      }
      void compute_hessian_sparsity(View<uno_int> row_indices, View<uno_int> column_indices, uno_int solver_indexing) const override {
         row_indices[0] = column_indices[0] = solver_indexing;
         row_indices[1] = column_indices[1] = 1 + solver_indexing;
      }
      void evaluate_lagrangian_hessian(const Vector<double>& /*x*/, double objective_multiplier, const Vector<double>& /*multipliers*/,
            View<double> hessian_values) const override {
         hessian_values[0] = 2. * objective_multiplier;
         hessian_values[1] = 2. * objective_multiplier;
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
      }
      void initial_dual_point(Vector<double>& multipliers) const override { multipliers.fill(0.); }
      void postprocess_solution(Iterate& /*iterate*/, Evaluations& /*evaluations*/) const override { }

      size_t number_jacobian_nonzeros() const override { return 2; }
      size_t number_hessian_nonzeros() const override { return 2; }

      size_t number_model_objective_evaluations() const override { return this->number_objective_evaluations; }
      size_t number_model_constraints_evaluations() const override { return this->number_constraint_evaluations; }
      size_t number_model_objective_gradient_evaluations() const override { return 0; }
      size_t number_model_jacobian_evaluations() const override { return 0; }
      size_t number_model_hessian_evaluations() const override { return 0; }
      void reset_number_evaluations() const override {
         this->number_objective_evaluations = this->number_constraint_evaluations = 0;
      }

      mutable size_t number_objective_evaluations{0};
      mutable size_t number_constraint_evaluations{0};

   private:
      const std::vector<uno_int> jacobian_rows{0, 0};
      const std::vector<uno_int> jacobian_columns{0, 1};
      const std::vector<double> variables_lower_bounds{0., -1.};
      const std::vector<double> variables_upper_bounds{INFINITY, 1.};
      const std::vector<double> constraints_lower_bounds{1.};
      const std::vector<double> constraints_upper_bounds{INFINITY};
      const std::vector<size_t> equality_constraints{};
      const std::vector<size_t> inequality_constraints{0};
      const IntegerRange linear_constraints{1};
      const IntegerRange nonlinear_constraints{0};
   };
} // namespace

#endif // UNO_TEST_BOUNDACTIVEMODEL_H
