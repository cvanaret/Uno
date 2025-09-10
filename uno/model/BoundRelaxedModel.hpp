// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BOUNDRELAXEDMODEL_H
#define UNO_BOUNDRELAXEDMODEL_H

#include "Model.hpp"

namespace uno {
   // forward declaration
   class Options;

   class BoundRelaxedModel: public Model {
   public:
      BoundRelaxedModel(const Model& original_model, const Options& options);

      // availability of linear operators
      [[nodiscard]] bool has_jacobian_operator() const override {
         return this->model.has_jacobian_operator();
      }

      [[nodiscard]] bool has_jacobian_transposed_operator() const override {
         return this->model.has_jacobian_transposed_operator();
      }

      [[nodiscard]] bool has_hessian_operator() const override {
         return this->model.has_hessian_operator();
      }

      [[nodiscard]] bool has_hessian_matrix() const override {
         return this->model.has_hessian_matrix();
      }

      // function evaluations
      [[nodiscard]] double evaluate_objective(const Vector<double>& x) const override {
         return this->model.evaluate_objective(x);
      }

      void evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const override {
         this->model.evaluate_constraints(x, constraints);
      }

      // dense objective gradient
      void evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const override {
         this->model.evaluate_objective_gradient(x, gradient);
      }

      // sparsity patterns of Jacobian and Hessian
      void compute_constraint_jacobian_sparsity(int* row_indices, int* column_indices, int solver_indexing,
            MatrixOrder matrix_order) const override {
         this->model.compute_constraint_jacobian_sparsity(row_indices, column_indices, solver_indexing, matrix_order);
      }

      void compute_hessian_sparsity(int* row_indices, int* column_indices, int solver_indexing) const override {
         this->model.compute_hessian_sparsity(row_indices, column_indices, solver_indexing);
      }

      void evaluate_constraint_jacobian(const Vector<double>& x, double* jacobian_values) const override {
         this->model.evaluate_constraint_jacobian(x, jacobian_values);
      }

      void evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
            double* hessian_values) const override {
         this->model.evaluate_lagrangian_hessian(x, objective_multiplier, multipliers, hessian_values);
      }

      void compute_jacobian_vector_product(const double* x, const double* vector, double* result) const override {
         this->model.compute_jacobian_vector_product(x, vector, result);
      }

      void compute_jacobian_transposed_vector_product(const double* x, const double* vector, double* result) const override {
         this->model.compute_jacobian_transposed_vector_product(x, vector, result);
      }

      void compute_hessian_vector_product(const double* x, const double* vector, double objective_multiplier,
            const Vector<double>& multipliers, double* result) const override {
         this->model.compute_hessian_vector_product(x, vector, objective_multiplier, multipliers, result);
      }

      // only these two functions are redefined
      [[nodiscard]] double variable_lower_bound(size_t variable_index) const override;
      [[nodiscard]] double variable_upper_bound(size_t variable_index) const override;
      [[nodiscard]] const SparseVector<size_t>& get_slacks() const override { return this->model.get_slacks(); }
      [[nodiscard]] const Vector<size_t>& get_fixed_variables() const override { return this->model.get_fixed_variables(); }

      [[nodiscard]] double constraint_lower_bound(size_t constraint_index) const override { return this->model.constraint_lower_bound(constraint_index); }
      [[nodiscard]] double constraint_upper_bound(size_t constraint_index) const override { return this->model.constraint_upper_bound(constraint_index); }
      [[nodiscard]] const Collection<size_t>& get_equality_constraints() const override { return this->model.get_equality_constraints(); }
      [[nodiscard]] const Collection<size_t>& get_inequality_constraints() const override { return this->model.get_inequality_constraints(); }
      [[nodiscard]] const Collection<size_t>& get_linear_constraints() const override { return this->model.get_linear_constraints(); }

      void initial_primal_point(Vector<double>& x) const override { this->model.initial_primal_point(x); }
      void initial_dual_point(Vector<double>& multipliers) const override { this->model.initial_dual_point(multipliers); }
      void postprocess_solution(Iterate& iterate) const override {
         this->model.postprocess_solution(iterate);
      }

      [[nodiscard]] size_t number_jacobian_nonzeros() const override { return this->model.number_jacobian_nonzeros(); }
      [[nodiscard]] size_t number_hessian_nonzeros() const override { return this->model.number_hessian_nonzeros(); }

   private:
      const Model& model;
      const double relaxation_factor;
   };
} // namespace

#endif // UNO_BOUNDRELAXEDMODEL_H