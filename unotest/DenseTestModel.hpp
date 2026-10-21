// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DENSETESTMODEL_H
#define UNO_DENSETESTMODEL_H

#include <functional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include "model/Model.hpp"
#include "linear_algebra/Vector.hpp"
#include "linear_algebra/View.hpp"
#include "symbolic/CollectionAdapter.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   // small dense NLP described by callbacks, for regression tests of the solver paths:
   //   min f(x) s.t. cl <= c(x) <= cu, xl <= x <= xu
   // The Jacobian is dense (row-major m x n) and the Lagrangian Hessian is the dense lower triangle
   // (row-major: (0,0), (1,0), (1,1), (2,0), ...)
   class DenseTestModel: public Model {
   public:
      using Objective = std::function<double(const Vector<double>&)>;
      using Gradient = std::function<void(const Vector<double>&, Vector<double>&)>;
      using Constraints = std::function<void(const Vector<double>&, Vector<double>&)>;
      using Jacobian = std::function<void(const Vector<double>&, double*)>;
      using Hessian = std::function<void(const Vector<double>&, double, const Vector<double>&, double*)>;

      DenseTestModel(std::string name, size_t number_variables, size_t number_constraints, ProblemType problem_type,
            Objective objective, Gradient gradient, Constraints constraints, Jacobian jacobian, Hessian hessian,
            std::vector<double> variables_lower_bounds, std::vector<double> variables_upper_bounds,
            std::vector<double> constraints_lower_bounds, std::vector<double> constraints_upper_bounds,
            std::vector<double> initial_point, std::vector<size_t> linear_constraints = {}):
               Model(std::move(name), number_variables, number_constraints, 1., -1., 0),
               problem_type(problem_type), objective(std::move(objective)), gradient(std::move(gradient)),
               constraints(std::move(constraints)), jacobian(std::move(jacobian)), hessian(std::move(hessian)),
               variables_lower_bounds(std::move(variables_lower_bounds)), variables_upper_bounds(std::move(variables_upper_bounds)),
               constraints_lower_bounds(std::move(constraints_lower_bounds)), constraints_upper_bounds(std::move(constraints_upper_bounds)),
               initial_point(std::move(initial_point)), linear_constraints(std::move(linear_constraints)) {
         for (size_t constraint_index = 0; constraint_index < number_constraints; ++constraint_index) {
            for (size_t variable_index = 0; variable_index < number_variables; ++variable_index) {
               this->jacobian_rows.push_back(static_cast<uno_int>(constraint_index));
               this->jacobian_columns.push_back(static_cast<uno_int>(variable_index));
            }
            bool is_linear = false;
            for (size_t linear_index: this->linear_constraints) {
               is_linear = is_linear || (linear_index == constraint_index);
            }
            if (!is_linear) {
               this->nonlinear_constraints.push_back(constraint_index);
            }
         }
         this->partition_constraints(this->equality_constraints, this->inequality_constraints);
      }
      // the constraint collections refer to members of the object: no copy or move (factories rely on guaranteed elision)
      DenseTestModel(const DenseTestModel&) = delete;
      DenseTestModel& operator=(const DenseTestModel&) = delete;

      [[nodiscard]] ProblemType get_problem_type() const override { return this->problem_type; }
      [[nodiscard]] bool has_jacobian_operator() const override { return false; }
      [[nodiscard]] bool has_jacobian_transposed_operator() const override { return false; }
      [[nodiscard]] bool has_hessian_operator() const override { return this->hessian_operator_available; }
      [[nodiscard]] bool has_hessian_matrix() const override { return true; }

      [[nodiscard]] double evaluate_objective(const Vector<double>& x) const override {
         ++this->number_objective_evaluations;
         return this->objective(x);
      }
      void evaluate_constraints(const Vector<double>& x, Vector<double>& c) const override {
         ++this->number_constraint_evaluations;
         this->constraints(x, c);
      }
      void evaluate_objective_gradient(const Vector<double>& x, Vector<double>& g) const override {
         ++this->number_gradient_evaluations;
         this->gradient(x, g);
      }
      [[nodiscard]] View<const uno_int> get_jacobian_row_indices() const override {
         return view(this->jacobian_rows.data(), this->jacobian_rows.size());
      }
      [[nodiscard]] View<const uno_int> get_jacobian_column_indices() const override {
         return view(this->jacobian_columns.data(), this->jacobian_columns.size());
      }
      void compute_hessian_sparsity(View<uno_int> row_indices, View<uno_int> column_indices, uno_int solver_indexing) const override {
         size_t nonzero_index = 0;
         for (size_t row = 0; row < this->number_variables; ++row) {
            for (size_t column = 0; column <= row; ++column) {
               row_indices[nonzero_index] = static_cast<uno_int>(row) + solver_indexing;
               column_indices[nonzero_index] = static_cast<uno_int>(column) + solver_indexing;
               ++nonzero_index;
            }
         }
      }
      void evaluate_jacobian(const Vector<double>& x, double* jacobian_values) const override {
         ++this->number_jacobian_evaluations;
         this->jacobian(x, jacobian_values);
      }
      void evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
            View<double> hessian_values) const override {
         ++this->number_hessian_evaluations;
         this->hessian(x, objective_multiplier, multipliers, hessian_values.data());
      }
      void compute_jacobian_vector_product(const double*, const double*, double*) const override {
         throw std::runtime_error("DenseTestModel: no Jacobian operator");
      }
      void compute_jacobian_transposed_vector_product(const double*, const double*, double*) const override {
         throw std::runtime_error("DenseTestModel: no Jacobian operator");
      }
      // Hessian-vector product computed from the dense lower triangle (available after enable_hessian_operator())
      void compute_hessian_vector_product(View<const double> x, View<const double> vector, double objective_multiplier,
            const Vector<double>& multipliers, View<double> result) const override {
         if (!this->hessian_operator_available) {
            throw std::runtime_error("DenseTestModel: no Hessian operator");
         }
         ++this->number_hessian_vector_products;
         Vector<double> primals(this->number_variables);
         for (size_t variable_index = 0; variable_index < this->number_variables; ++variable_index) {
            primals[variable_index] = x[variable_index];
         }
         std::vector<double> lower_triangle(this->number_hessian_nonzeros());
         this->hessian(primals, objective_multiplier, multipliers, lower_triangle.data());
         for (size_t variable_index = 0; variable_index < this->number_variables; ++variable_index) {
            result[variable_index] = 0.;
         }
         size_t nonzero_index = 0;
         for (size_t row = 0; row < this->number_variables; ++row) {
            for (size_t column = 0; column <= row; ++column) {
               result[row] += lower_triangle[nonzero_index] * vector[column];
               if (row != column) {
                  result[column] += lower_triangle[nonzero_index] * vector[row];
               }
               ++nonzero_index;
            }
         }
      }

      // expose a Hessian operator in addition to the Hessian matrix (like the AMPL interface)
      void enable_hessian_operator() { this->hessian_operator_available = true; }
      [[nodiscard]] size_t number_model_hessian_vector_products() const { return this->number_hessian_vector_products; }

      [[nodiscard]] const std::vector<double>& get_variables_lower_bounds() const override { return this->variables_lower_bounds; }
      [[nodiscard]] const std::vector<double>& get_variables_upper_bounds() const override { return this->variables_upper_bounds; }
      [[nodiscard]] const std::vector<double>& get_constraints_lower_bounds() const override { return this->constraints_lower_bounds; }
      [[nodiscard]] const std::vector<double>& get_constraints_upper_bounds() const override { return this->constraints_upper_bounds; }
      [[nodiscard]] const std::vector<size_t>& get_equality_constraints() const override { return this->equality_constraints; }
      [[nodiscard]] const std::vector<size_t>& get_inequality_constraints() const override { return this->inequality_constraints; }
      [[nodiscard]] const Collection<size_t>& get_linear_constraints() const override { return this->linear_constraints_collection; }
      [[nodiscard]] const Collection<size_t>& get_nonlinear_constraints() const override { return this->nonlinear_constraints_collection; }

      void initial_primal_point(Vector<double>& x) const override {
         for (size_t variable_index = 0; variable_index < this->number_variables; ++variable_index) {
            x[variable_index] = this->initial_point[variable_index];
         }
      }
      void initial_dual_point(Vector<double>& multipliers) const override { multipliers.fill(0.); }
      void postprocess_solution(Iterate& /*iterate*/, Evaluations& /*evaluations*/) const override { }

      [[nodiscard]] size_t number_jacobian_nonzeros() const override { return this->jacobian_rows.size(); }
      [[nodiscard]] size_t number_hessian_nonzeros() const override { return this->number_variables * (this->number_variables + 1) / 2; }

      [[nodiscard]] size_t number_model_objective_evaluations() const override { return this->number_objective_evaluations; }
      [[nodiscard]] size_t number_model_constraints_evaluations() const override { return this->number_constraint_evaluations; }
      [[nodiscard]] size_t number_model_objective_gradient_evaluations() const override { return this->number_gradient_evaluations; }
      [[nodiscard]] size_t number_model_jacobian_evaluations() const override { return this->number_jacobian_evaluations; }
      [[nodiscard]] size_t number_model_hessian_evaluations() const override { return this->number_hessian_evaluations; }
      void reset_number_evaluations() const override {
         this->number_objective_evaluations = this->number_constraint_evaluations = this->number_gradient_evaluations =
            this->number_jacobian_evaluations = this->number_hessian_evaluations = 0;
      }

   private:
      const ProblemType problem_type;
      const Objective objective;
      const Gradient gradient;
      const Constraints constraints;
      const Jacobian jacobian;
      const Hessian hessian;
      const std::vector<double> variables_lower_bounds, variables_upper_bounds;
      const std::vector<double> constraints_lower_bounds, constraints_upper_bounds;
      const std::vector<double> initial_point;
      const std::vector<size_t> linear_constraints;
      std::vector<size_t> nonlinear_constraints{};
      std::vector<size_t> equality_constraints{}, inequality_constraints{};
      std::vector<uno_int> jacobian_rows{}, jacobian_columns{};
      const CollectionAdapter<std::vector<size_t>> linear_constraints_collection{this->linear_constraints};
      const CollectionAdapter<std::vector<size_t>> nonlinear_constraints_collection{this->nonlinear_constraints};
      mutable size_t number_objective_evaluations{0};
      mutable size_t number_constraint_evaluations{0};
      mutable size_t number_gradient_evaluations{0};
      mutable size_t number_jacobian_evaluations{0};
      mutable size_t number_hessian_evaluations{0};
      bool hessian_operator_available{false};
      mutable size_t number_hessian_vector_products{0};
   };
} // namespace

#endif // UNO_DENSETESTMODEL_H