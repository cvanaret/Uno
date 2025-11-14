// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRIMALDUALINTERIORPOINTPROBLEM_H
#define UNO_PRIMALDUALINTERIORPOINTPROBLEM_H

#include "../BarrierProblem.hpp"
#include "../InteriorPointParameters.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "symbolic/Range.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   class PrimalDualInteriorPointProblem : public BarrierProblem {
   public:
      PrimalDualInteriorPointProblem(const OptimizationProblem& problem, const InteriorPointParameters &parameters);

      [[nodiscard]] double get_objective_multiplier() const override;
      void set_barrier_parameter(double barrier_parameter) override;
      void generate_initial_iterate(Iterate& initial_iterate) const override;

      // constraint evaluations
      void evaluate_constraints(Iterate& iterate, Vector<double>& constraints) const override;

      // dense objective gradient
      void evaluate_objective_gradient(Iterate& iterate, double* objective_gradient) const override;

      // sparsity patterns of Jacobian and Hessian
      void compute_constraint_jacobian_sparsity(uno_int* row_indices, uno_int* column_indices, uno_int solver_indexing,
         MatrixOrder matrix_order) const override;
      void compute_hessian_sparsity(const HessianModel& hessian_model, uno_int* row_indices,
         uno_int* column_indices, uno_int solver_indexing) const override;

      // numerical evaluations of Jacobian and Hessian
      [[nodiscard]] size_t number_jacobian_nonzeros() const override;
      [[nodiscard]] bool has_curvature(const HessianModel& hessian_model) const override;
      [[nodiscard]] size_t number_hessian_nonzeros(const HessianModel& hessian_model) const override;
      void evaluate_constraint_jacobian(Iterate& iterate, double* jacobian_values) const override;
      void evaluate_lagrangian_gradient(LagrangianGradient<double>& lagrangian_gradient,
         const EvaluationSpace& evaluation_space, Iterate& iterate) const override;
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

      void assemble_primal_dual_direction(const Iterate& current_iterate, const Vector<double>& solution, Direction& direction) const override;

      [[nodiscard]] double push_variable_to_interior(double variable_value, double lower_bound, double upper_bound) const;
      [[nodiscard]] double dual_regularization_factor() const override;
      [[nodiscard]] double compute_barrier_term_directional_derivative(const Iterate& current_iterate,
         const Vector<double>& primal_direction) const;
      void postprocess_iterate(Iterate& iterate) const override;
      [[nodiscard]] double compute_centrality_error(const Vector<double>& primals, const Multipliers& multipliers,
         double shift) const override;

      // progress measures
      void set_auxiliary_measure(Iterate& iterate) const override;
      [[nodiscard]] double compute_predicted_auxiliary_reduction(const Iterate& current_iterate,
         const Vector<double>& primal_direction, double step_length) const override;

   protected:
      const OptimizationProblem& first_reformulation;
      double barrier_parameter{INF<double>};
      const InteriorPointParameters& parameters;
      const Vector<size_t> fixed_variables{};
      const ForwardRange equality_constraints;
      const ForwardRange inequality_constraints{0};

      void compute_bound_dual_direction(const Iterate& current_iterate, Direction& direction) const;
      [[nodiscard]] double primal_fraction_to_boundary(const Vector<double>& current_primals, const Vector<double>& primal_direction,
         double tau) const;
      [[nodiscard]] double dual_fraction_to_boundary(const Multipliers& current_multipliers, const Multipliers& direction_multipliers,
         double tau) const;
   };
} // namespace

#endif // UNO_PRIMALDUALINTERIORPOINTPROBLEM_H