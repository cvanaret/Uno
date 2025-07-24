// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SCALEDMODEL_H
#define UNO_SCALEDMODEL_H

#include <memory>
#include "Model.hpp"
#include "linear_algebra/Vector.hpp"
#include "preprocessing/Scaling.hpp"

namespace uno {
   // forward declaration
   class Options;

   class ScaledModel: public Model {
   public:
      ScaledModel(std::unique_ptr<Model> original_model, Iterate& initial_iterate, const Options& options);

      // function evaluations
      [[nodiscard]] double evaluate_objective(const Vector<double>& x) const override;
      void evaluate_constraints(const Vector<double>& x, std::vector<double>& constraints) const override;

      // dense objective gradient
      void evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const override;

      // structures of Jacobian and Hessian
      void compute_jacobian_structure(Vector<size_t>& row_indices, Vector<size_t>& column_indices) const override;
      void compute_hessian_structure(Vector<size_t>& row_indices, Vector<size_t>& column_indices) const override;

      // numerical evaluations of Jacobian and Hessian
      void evaluate_constraint_jacobian(const Vector<double>& x, Vector<double>& jacobian_values) const override;
      void evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         Vector<double>& hessian_values) const override;
      void compute_hessian_vector_product(const double* vector, double objective_multiplier, const Vector<double>& multipliers,
         double* result) const override;

      [[nodiscard]] double variable_lower_bound(size_t variable_index) const override;
      [[nodiscard]] double variable_upper_bound(size_t variable_index) const override;
      [[nodiscard]] const Collection<size_t>& get_lower_bounded_variables() const override;
      [[nodiscard]] const Collection<size_t>& get_upper_bounded_variables() const override;
      [[nodiscard]] const SparseVector<size_t>& get_slacks() const override;
      [[nodiscard]] const Collection<size_t>& get_single_lower_bounded_variables() const override;
      [[nodiscard]] const Collection<size_t>& get_single_upper_bounded_variables() const override;
      [[nodiscard]] const Vector<size_t>& get_fixed_variables() const override;

      [[nodiscard]] double constraint_lower_bound(size_t constraint_index) const override;
      [[nodiscard]] double constraint_upper_bound(size_t constraint_index) const override;
      [[nodiscard]] const Collection<size_t>& get_equality_constraints() const override;
      [[nodiscard]] const Collection<size_t>& get_inequality_constraints() const override;
      [[nodiscard]] const Collection<size_t>& get_linear_constraints() const override;

      void initial_primal_point(Vector<double>& x) const override;
      void initial_dual_point(Vector<double>& multipliers) const override;
      void postprocess_solution(Iterate& iterate, IterateStatus termination_status) const override;

      [[nodiscard]] size_t number_jacobian_nonzeros() const override;
      [[nodiscard]] size_t number_hessian_nonzeros() const override;

   private:
      const std::unique_ptr<Model> model{};
      Scaling scaling;
      mutable Vector<double> scaled_multipliers{};
   };
} // namespace

#endif // UNO_SCALEDMODEL_H