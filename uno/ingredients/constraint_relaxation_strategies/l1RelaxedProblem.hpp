// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_L1RELAXEDPROBLEM_H
#define UNO_L1RELAXEDPROBLEM_H

#include <functional>
#include "optimization/OptimizationProblem.hpp"
#include "symbolic/Concatenation.hpp"

namespace uno {
   class l1RelaxedProblem: public OptimizationProblem {
   public:
      // constructor with proximal term
      l1RelaxedProblem(const Model& model, double objective_multiplier, double constraint_violation_coefficient, double proximal_coefficient,
         double const* proximal_center);
      // constructor without proximal term
      l1RelaxedProblem(const Model& model, double objective_multiplier, double constraint_violation_coefficient);

      [[nodiscard]] double get_objective_multiplier() const override;
      void evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const override;
      void evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const override;
      void evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const override;
      void evaluate_lagrangian_hessian(Statistics& statistics, HessianModel& hessian_model, const Vector<double>& primal_variables,
         const Multipliers& multipliers, SymmetricMatrix<size_t, double>& hessian) const override;
      void compute_hessian_vector_product(HessianModel& hessian_model, const Vector<double>& vector, const Multipliers& multipliers,
         Vector<double>& result) const override;

      [[nodiscard]] double variable_lower_bound(size_t variable_index) const override;
      [[nodiscard]] double variable_upper_bound(size_t variable_index) const override;
      [[nodiscard]] const Collection<size_t>& get_lower_bounded_variables() const override;
      [[nodiscard]] const Collection<size_t>& get_upper_bounded_variables() const override;
      [[nodiscard]] const Collection<size_t>& get_single_lower_bounded_variables() const override;
      [[nodiscard]] const Collection<size_t>& get_single_upper_bounded_variables() const override;
      [[nodiscard]] const Vector<size_t>& get_fixed_variables() const override;

      [[nodiscard]] double constraint_lower_bound(size_t constraint_index) const override;
      [[nodiscard]] double constraint_upper_bound(size_t constraint_index) const override;
      [[nodiscard]] const Collection<size_t>& get_equality_constraints() const override;
      [[nodiscard]] const Collection<size_t>& get_inequality_constraints() const override;

      [[nodiscard]] size_t number_objective_gradient_nonzeros() const override;
      [[nodiscard]] size_t number_jacobian_nonzeros() const override;
      [[nodiscard]] size_t number_hessian_nonzeros(const HessianModel& hessian_model) const override;

      void evaluate_lagrangian_gradient(LagrangianGradient<double>& lagrangian_gradient, Iterate& iterate, const Multipliers& multipliers) const override;
      [[nodiscard]] double complementarity_error(const Vector<double>& primals, const std::vector<double>& constraints,
         const Multipliers& multipliers, double shift_value, Norm residual_norm) const override;

      // parameterization
      void set_proximal_multiplier(double new_proximal_coefficient);
      void set_proximal_center(double const* new_proximal_center);
      void set_elastic_variable_values(Iterate& iterate, const std::function<void(Iterate&, size_t, size_t, double)>& elastic_setting_function) const;

   protected:
      const size_t number_elastic_variables;
      const double objective_multiplier;
      const double constraint_violation_coefficient;
      double proximal_coefficient;
      double const* proximal_center;
      const Concatenation<const Collection<size_t>&, ForwardRange> lower_bounded_variables; // model variables + elastic variables
      const Concatenation<const Collection<size_t>&, ForwardRange> single_lower_bounded_variables; // model variables + elastic variables
   };
} // namespace

#endif // UNO_L1RELAXEDPROBLEM_H
