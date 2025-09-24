// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_EXPONENTIALBARRIERPROBLEM_H
#define UNO_EXPONENTIALBARRIERPROBLEM_H

#include "../InteriorPointParameters.hpp"
#include "../BarrierProblem.hpp"
#include "symbolic/Range.hpp"

namespace uno {
   class ExponentialBarrierProblem : public BarrierProblem {
   public:
      ExponentialBarrierProblem(const OptimizationProblem& problem, double barrier_parameter,
         const InteriorPointParameters &parameters);

      // function evaluations
      [[nodiscard]] double get_objective_multiplier() const override;

      // constraint evaluations
      void evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const override;

      // dense objective gradient
      void evaluate_objective_gradient(Iterate& iterate, const EvaluationSpace& evaluation_space, double* objective_gradient) const override;

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
         const InequalityHandlingMethod& inequality_handling_method, const EvaluationSpace& evaluation_space, Iterate& iterate) const override;
      void evaluate_lagrangian_hessian(Statistics& statistics, HessianModel& hessian_model, const Vector<double>& primal_variables,
         const Multipliers& multipliers, double* hessian_values) const override;
      void compute_hessian_vector_product(HessianModel& hessian_model, const double* vector,
         const Multipliers& multipliers, double* result) const override;

      [[nodiscard]] double variable_lower_bound(size_t variable_index) const override;
      [[nodiscard]] double variable_upper_bound(size_t variable_index) const override;
      [[nodiscard]] const Collection<size_t>& get_lower_bounded_variables() const override;
      [[nodiscard]] const Collection<size_t>& get_upper_bounded_variables() const override;
      [[nodiscard]] const Collection<size_t>& get_single_lower_bounded_variables() const override;
      [[nodiscard]] const Collection<size_t>& get_single_upper_bounded_variables() const override;
      [[nodiscard]] const Vector<size_t>& get_fixed_variables() const override;
      [[nodiscard]] const Collection<size_t>& get_primal_regularization_variables() const override;

      [[nodiscard]] double constraint_lower_bound(size_t constraint_index) const override;
      [[nodiscard]] double constraint_upper_bound(size_t constraint_index) const override;
      [[nodiscard]] const Collection<size_t>& get_equality_constraints() const override;
      [[nodiscard]] const Collection<size_t>& get_inequality_constraints() const override;
      [[nodiscard]] const Collection<size_t>& get_dual_regularization_constraints() const override;

      void assemble_primal_dual_direction(const Iterate& current_iterate, const Vector<double>& solution,
         Direction& direction) const override;

      [[nodiscard]] double push_variable_to_interior(double variable_value, double lower_bound, double upper_bound) const;
      void set_auxiliary_measure(Iterate& iterate) const;
      [[nodiscard]] double complementarity_error(const Vector<double>& primals, const std::vector<double>& constraints,
         const Multipliers& multipliers, double shift_value, Norm residual_norm) const override;
      [[nodiscard]] double dual_regularization_factor() const override;
      [[nodiscard]] double compute_barrier_term_directional_derivative(const Iterate& current_iterate,
         const Vector<double>& primal_direction) const;
      void postprocess_iterate(Iterate& iterate) const;
      [[nodiscard]] double compute_centrality_error(const Vector<double>& primals, const Multipliers& multipliers,
         double barrier_parameter) const override;

   protected:
      const OptimizationProblem& reformulated_problem;
      const double barrier_parameter;
      const InteriorPointParameters& parameters;
      const Vector<size_t> fixed_variables{};
      const ForwardRange inequality_constraints{0};

      [[nodiscard]] static size_t count_number_variables(const OptimizationProblem& problem);
   };
} // namespace

#endif // UNO_EXPONENTIALBARRIERPROBLEM_H