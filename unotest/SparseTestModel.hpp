// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SPARSETESTMODEL_H
#define UNO_SPARSETESTMODEL_H

#include <array>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include "model/Model.hpp"
#include "linear_algebra/Vector.hpp"
#include "linear_algebra/View.hpp"
#include "symbolic/CollectionAdapter.hpp"

namespace uno {
   // small NLP min f(x) s.t. cl <= c(x) <= cu, xl <= x <= xu with sparse (0-based, COO) Jacobian and Hessian (lower
   // triangle), described by callbacks
   struct SparseTestProblem {
      std::string name;
      size_t number_variables;
      size_t number_constraints;
      std::function<double(const Vector<double>&)> objective;
      std::function<void(const Vector<double>&, Vector<double>&)> gradient;
      std::function<void(const Vector<double>&, Vector<double>&)> constraints;
      std::vector<uno_int> jacobian_rows, jacobian_columns;
      std::function<void(const Vector<double>&, double*)> jacobian;
      std::vector<uno_int> hessian_rows, hessian_columns;
      // Hessian of the Lagrangian ρ f(x) - y^T c(x)
      std::function<void(const Vector<double>&, double, const Vector<double>&, double*)> hessian;
      std::vector<double> variables_lower_bounds, variables_upper_bounds;
      std::vector<double> constraints_lower_bounds, constraints_upper_bounds;
      std::vector<double> initial_point;
      std::vector<size_t> linear_constraints;
   };

   class SparseTestModel: public Model {
   public:
      explicit SparseTestModel(SparseTestProblem problem):
            Model(problem.name, problem.number_variables, problem.number_constraints, 1., -1., 0),
            problem(std::move(problem)) {
         std::vector<bool> is_linear(this->number_constraints, false);
         for (size_t constraint_index: this->problem.linear_constraints) {
            is_linear[constraint_index] = true;
         }
         for (size_t constraint_index = 0; constraint_index < this->number_constraints; ++constraint_index) {
            if (!is_linear[constraint_index]) {
               this->nonlinear_constraints.push_back(constraint_index);
            }
         }
         this->partition_constraints(this->equality_constraints, this->inequality_constraints);
      }
      // the constraint collections refer to members of the object
      SparseTestModel(const SparseTestModel&) = delete;
      SparseTestModel& operator=(const SparseTestModel&) = delete;

      [[nodiscard]] ProblemType get_problem_type() const override { return ProblemType::NONLINEAR; }
      [[nodiscard]] bool has_jacobian_operator() const override { return false; }
      [[nodiscard]] bool has_jacobian_transposed_operator() const override { return false; }
      [[nodiscard]] bool has_hessian_operator() const override { return false; }
      [[nodiscard]] bool has_hessian_matrix() const override { return true; }

      [[nodiscard]] double evaluate_objective(const Vector<double>& x) const override {
         ++this->number_evaluations[0];
         return this->problem.objective(x);
      }
      void evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const override {
         ++this->number_evaluations[1];
         this->problem.constraints(x, constraints);
      }
      void evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const override {
         ++this->number_evaluations[2];
         this->problem.gradient(x, gradient);
      }
      [[nodiscard]] View<const uno_int> get_jacobian_row_indices() const override {
         return view(this->problem.jacobian_rows.data(), this->problem.jacobian_rows.size());
      }
      [[nodiscard]] View<const uno_int> get_jacobian_column_indices() const override {
         return view(this->problem.jacobian_columns.data(), this->problem.jacobian_columns.size());
      }
      void compute_hessian_sparsity(View<uno_int> row_indices, View<uno_int> column_indices, uno_int solver_indexing) const override {
         for (size_t nonzero_index = 0; nonzero_index < this->problem.hessian_rows.size(); ++nonzero_index) {
            row_indices[nonzero_index] = this->problem.hessian_rows[nonzero_index] + solver_indexing;
            column_indices[nonzero_index] = this->problem.hessian_columns[nonzero_index] + solver_indexing;
         }
      }
      void evaluate_jacobian(const Vector<double>& x, double* jacobian_values) const override {
         ++this->number_evaluations[3];
         this->problem.jacobian(x, jacobian_values);
      }
      void evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
            View<double> hessian_values) const override {
         ++this->number_evaluations[4];
         this->problem.hessian(x, objective_multiplier, multipliers, hessian_values.data());
      }
      void compute_jacobian_vector_product(const double*, const double*, double*) const override {
         throw std::runtime_error("SparseTestModel: no Jacobian operator");
      }
      void compute_jacobian_transposed_vector_product(const double*, const double*, double*) const override {
         throw std::runtime_error("SparseTestModel: no Jacobian operator");
      }
      void compute_hessian_vector_product(View<const double>, View<const double>, double, const Vector<double>&,
            View<double>) const override {
         throw std::runtime_error("SparseTestModel: no Hessian operator");
      }

      [[nodiscard]] const std::vector<double>& get_variables_lower_bounds() const override { return this->problem.variables_lower_bounds; }
      [[nodiscard]] const std::vector<double>& get_variables_upper_bounds() const override { return this->problem.variables_upper_bounds; }
      [[nodiscard]] const std::vector<double>& get_constraints_lower_bounds() const override { return this->problem.constraints_lower_bounds; }
      [[nodiscard]] const std::vector<double>& get_constraints_upper_bounds() const override { return this->problem.constraints_upper_bounds; }
      [[nodiscard]] const std::vector<size_t>& get_equality_constraints() const override { return this->equality_constraints; }
      [[nodiscard]] const std::vector<size_t>& get_inequality_constraints() const override { return this->inequality_constraints; }
      [[nodiscard]] const Collection<size_t>& get_linear_constraints() const override { return this->linear_constraints_collection; }
      [[nodiscard]] const Collection<size_t>& get_nonlinear_constraints() const override { return this->nonlinear_constraints_collection; }

      void initial_primal_point(Vector<double>& x) const override {
         for (size_t variable_index = 0; variable_index < this->number_variables; ++variable_index) {
            x[variable_index] = this->problem.initial_point[variable_index];
         }
      }
      void initial_dual_point(Vector<double>& multipliers) const override { multipliers.fill(0.); }
      void postprocess_solution(Iterate& /*iterate*/, Evaluations& /*evaluations*/) const override { }

      [[nodiscard]] size_t number_jacobian_nonzeros() const override { return this->problem.jacobian_rows.size(); }
      [[nodiscard]] size_t number_hessian_nonzeros() const override { return this->problem.hessian_rows.size(); }

      [[nodiscard]] size_t number_model_objective_evaluations() const override { return this->number_evaluations[0]; }
      [[nodiscard]] size_t number_model_constraints_evaluations() const override { return this->number_evaluations[1]; }
      [[nodiscard]] size_t number_model_objective_gradient_evaluations() const override { return this->number_evaluations[2]; }
      [[nodiscard]] size_t number_model_jacobian_evaluations() const override { return this->number_evaluations[3]; }
      [[nodiscard]] size_t number_model_hessian_evaluations() const override { return this->number_evaluations[4]; }
      void reset_number_evaluations() const override { this->number_evaluations = {}; }

   private:
      const SparseTestProblem problem;
      std::vector<size_t> nonlinear_constraints{};
      std::vector<size_t> equality_constraints{}, inequality_constraints{};
      const CollectionAdapter<std::vector<size_t>> linear_constraints_collection{this->problem.linear_constraints};
      const CollectionAdapter<std::vector<size_t>> nonlinear_constraints_collection{this->nonlinear_constraints};
      mutable std::array<size_t, 5> number_evaluations{};
   };
} // namespace

#endif // UNO_SPARSETESTMODEL_H
