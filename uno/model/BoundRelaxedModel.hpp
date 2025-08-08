// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BOUNDRELAXEDMODEL_H
#define UNO_BOUNDRELAXEDMODEL_H

#include <memory>
#include "Model.hpp"

namespace uno {
   // forward declaration
   class Options;

   class BoundRelaxedModel: public Model {
   public:
      BoundRelaxedModel(std::unique_ptr<Model> original_model, const Options& options);

      // function evaluations
      [[nodiscard]] double evaluate_objective(const Vector<double>& x) const override { return this->model->evaluate_objective(x); }
      void evaluate_constraints(const Vector<double>& x, std::vector<double>& constraints) const override {
         this->model->evaluate_constraints(x, constraints);
      }

      // dense objective gradient
      void evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const override {
         this->model->evaluate_objective_gradient(x, gradient);
      }

      // sparsity patterns of Jacobian and Hessian
      void compute_constraint_jacobian_sparsity(size_t* row_indices, size_t* column_indices, size_t solver_indexing,
            MatrixOrder matrix_order) const override {
         this->model->compute_constraint_jacobian_sparsity(row_indices, column_indices, solver_indexing, matrix_order);
      }

      void compute_hessian_sparsity(size_t* row_indices, size_t* column_indices, size_t solver_indexing) const override {
         this->model->compute_hessian_sparsity(row_indices, column_indices, solver_indexing);
      }

      void evaluate_constraint_jacobian(const Vector<double>& x, double* jacobian_values) const override {
         this->model->evaluate_constraint_jacobian(x, jacobian_values);
      }

      void evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
            double* hessian_values) const override {
         this->model->evaluate_lagrangian_hessian(x, objective_multiplier, multipliers, hessian_values);
      }

      void compute_hessian_vector_product(const double* vector, double objective_multiplier, const Vector<double>& multipliers,
            double* result) const override {
         this->model->compute_hessian_vector_product(vector, objective_multiplier, multipliers, result);
      }

      // only these two functions are redefined
      [[nodiscard]] double variable_lower_bound(size_t variable_index) const override;
      [[nodiscard]] double variable_upper_bound(size_t variable_index) const override;
      [[nodiscard]] const Collection<size_t>& get_lower_bounded_variables() const override { return this->model->get_lower_bounded_variables(); }
      [[nodiscard]] const Collection<size_t>& get_upper_bounded_variables() const override { return this->model->get_upper_bounded_variables(); }
      [[nodiscard]] const SparseVector<size_t>& get_slacks() const override { return this->model->get_slacks(); }
      [[nodiscard]] const Collection<size_t>& get_single_lower_bounded_variables() const override { return this->model->get_single_lower_bounded_variables(); }
      [[nodiscard]] const Collection<size_t>& get_single_upper_bounded_variables() const override { return this->model->get_single_upper_bounded_variables(); }
      [[nodiscard]] const Vector<size_t>& get_fixed_variables() const override { return this->model->get_fixed_variables(); }

      [[nodiscard]] double constraint_lower_bound(size_t constraint_index) const override { return this->model->constraint_lower_bound(constraint_index); }
      [[nodiscard]] double constraint_upper_bound(size_t constraint_index) const override { return this->model->constraint_upper_bound(constraint_index); }
      [[nodiscard]] const Collection<size_t>& get_equality_constraints() const override { return this->model->get_equality_constraints(); }
      [[nodiscard]] const Collection<size_t>& get_inequality_constraints() const override { return this->model->get_inequality_constraints(); }
      [[nodiscard]] const Collection<size_t>& get_linear_constraints() const override { return this->model->get_linear_constraints(); }

      void initial_primal_point(Vector<double>& x) const override { this->model->initial_primal_point(x); }
      void initial_dual_point(Vector<double>& multipliers) const override { this->model->initial_dual_point(multipliers); }
      void postprocess_solution(Iterate& iterate, IterateStatus termination_status) const override {
         this->model->postprocess_solution(iterate, termination_status);
      }

      [[nodiscard]] size_t number_jacobian_nonzeros() const override { return this->model->number_jacobian_nonzeros(); }
      [[nodiscard]] size_t number_hessian_nonzeros() const override { return this->model->number_hessian_nonzeros(); }

   private:
      const std::unique_ptr<Model> model{};
      const double relaxation_factor;
   };
} // namespace

#endif // UNO_BOUNDRELAXEDMODEL_H