// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_L1RELAXEDPROBLEM_H
#define UNO_L1RELAXEDPROBLEM_H

#include <functional>
#include "optimization/OptimizationProblem.hpp"

namespace uno {
   class l1RelaxedProblem: public OptimizationProblem {
   public:
      l1RelaxedProblem(const Model& model, double objective_multiplier, double constraint_violation_coefficient,
         double proximal_coefficient, double const* proximal_center);
      ~l1RelaxedProblem() override = default;

      [[nodiscard]] double get_objective_multiplier() const override;

      // constraint evaluations
      void evaluate_constraints(Iterate& iterate, Vector<double>& constraints) const override;

      // dense objective gradient
      void evaluate_objective_gradient(Iterate& iterate, double* objective_gradient) const override;

      // sparsity patterns of Jacobian and Hessian
      [[nodiscard]] size_t number_jacobian_nonzeros() const override;
      [[nodiscard]] bool has_curvature(const HessianModel& hessian_model) const override;
      [[nodiscard]] size_t number_hessian_nonzeros(const HessianModel& hessian_model) const override;
      void compute_constraint_jacobian_sparsity(int* row_indices, int* column_indices, int solver_indexing,
         MatrixOrder matrix_order) const override;
      void compute_hessian_sparsity(const HessianModel& hessian_model, int* row_indices,
         int* column_indices, int solver_indexing) const override;

      // numerical evaluations of Jacobian and Hessian
      void evaluate_constraint_jacobian(Iterate& iterate, double* jacobian_values) const override;
      void evaluate_lagrangian_gradient(LagrangianGradient<double>& lagrangian_gradient,
         const InequalityHandlingMethod& inequality_handling_method, Iterate& iterate) const override;
      void evaluate_lagrangian_hessian(Statistics& statistics, HessianModel& hessian_model, const Vector<double>& primal_variables,
         const Multipliers& multipliers, double* hessian_values) const override;
      void compute_hessian_vector_product(HessianModel& hessian_model, const double* x, const double* vector,
         const Multipliers& multipliers, double* result) const override;

      [[nodiscard]] double variable_lower_bound(size_t variable_index) const override;
      [[nodiscard]] double variable_upper_bound(size_t variable_index) const override;
      [[nodiscard]] const Vector<size_t>& get_fixed_variables() const override;
      // [[nodiscard]] virtual const Collection<size_t>& get_primal_regularization_variables() const;

      [[nodiscard]] double constraint_lower_bound(size_t constraint_index) const override;
      [[nodiscard]] double constraint_upper_bound(size_t constraint_index) const override;
      [[nodiscard]] const Collection<size_t>& get_equality_constraints() const override;
      [[nodiscard]] const Collection<size_t>& get_inequality_constraints() const override;
      [[nodiscard]] const Collection<size_t>& get_dual_regularization_constraints() const override;

      [[nodiscard]] SolutionStatus check_first_order_convergence(const Iterate& current_iterate, double primal_tolerance,
         double dual_tolerance) const;

      void set_elastic_variable_values(Iterate& iterate, const std::function<void(Iterate&, size_t, size_t, double)>& elastic_setting_function) const;

   protected:
      const size_t number_elastic_variables;
      const double objective_multiplier;
      const double constraint_violation_coefficient;
      const double proximal_coefficient;
      double const* proximal_center;
      const ForwardRange dual_regularization_constraints{0};
   };
} // namespace

#endif // UNO_L1RELAXEDPROBLEM_H