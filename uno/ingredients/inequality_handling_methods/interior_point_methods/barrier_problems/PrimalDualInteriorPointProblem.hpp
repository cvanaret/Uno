// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRIMALDUALINTERIORPOINTPROBLEM_H
#define UNO_PRIMALDUALINTERIORPOINTPROBLEM_H

#include "optimization/OptimizationProblem.hpp"
#include "linear_algebra/Vector.hpp"
#include "symbolic/Range.hpp"

namespace uno {
   // forward declarations
   struct InteriorPointParameters;
   class Parameterization;

   class PrimalDualInteriorPointProblem : public OptimizationProblem {
   public:
      PrimalDualInteriorPointProblem(const OptimizationProblem& problem, const InteriorPointParameters& parameters,
         const Parameterization& parameterization);
      std::unique_ptr<OptimizationProblem> clone() const override;

      [[nodiscard]] double get_objective_multiplier() const override;
      void generate_initial_iterate(Iterate& initial_iterate, Evaluations& evaluations) const override;

      // sparsity patterns of Jacobian and Hessian
      [[nodiscard]] size_t number_jacobian_nonzeros() const override;
      [[nodiscard]] bool has_curvature(const HessianModel& hessian_model) const override;
      [[nodiscard]] size_t number_hessian_nonzeros(const HessianModel& hessian_model) const override;
      void compute_jacobian_sparsity(uno_int* row_indices, uno_int* column_indices, uno_int row_offset, uno_int column_offset,
         uno_int solver_indexing, MatrixOrder matrix_order) const override;
      void compute_hessian_sparsity(const HessianModel& hessian_model, uno_int* row_indices,
         uno_int* column_indices, uno_int solver_indexing) const override;

      // numerical evaluations of constraints, objective gradient, Jacobian and Hessian
      void evaluate_constraints(const Iterate& iterate, double* constraints, Evaluations& evaluations) const override;
      void evaluate_objective_gradient(const Iterate& iterate, double* objective_gradient, Evaluations& evaluations) const override;
      void evaluate_jacobian(const Vector<double>& primals, double* jacobian_values, Evaluations& evaluations) const override;
      void evaluate_lagrangian_gradient(const Iterate& iterate, Evaluations& evaluations,
         Vector<double>& lagrangian_gradient) const override;
      void evaluate_lagrangian_hessian(Statistics& statistics, HessianModel& hessian_model, const Vector<double>& primal_variables,
         const Multipliers& multipliers, double* hessian_values) const override;

      // linear operators
      void compute_jacobian_vector_product(const double* vector, double* result, const Evaluations& evaluations) const override;
      void compute_jacobian_transposed_vector_product(const double* vector, double* result, const Evaluations& evaluations) const override;
      void compute_hessian_vector_product(HessianModel& hessian_model, const double* x, const double* vector,
         const Multipliers& multipliers, double* result) const override;

      [[nodiscard]] const std::vector<double>& get_variables_lower_bounds() const override;
      [[nodiscard]] const std::vector<double>& get_variables_upper_bounds() const override;
      [[nodiscard]] const Vector<size_t>& get_fixed_variables() const override;
      // [[nodiscard]] virtual const Collection<size_t>& get_primal_regularization_variables() const;

      [[nodiscard]] const std::vector<double>& get_constraints_lower_bounds() const override;
      [[nodiscard]] const std::vector<double>& get_constraints_upper_bounds() const override;
      [[nodiscard]] const Collection<size_t>& get_equality_constraints() const override;
      [[nodiscard]] const Collection<size_t>& get_inequality_constraints() const override;
      [[nodiscard]] const Collection<size_t>& get_dual_regularization_constraints() const override;

      void assemble_primal_dual_direction(const Iterate& current_iterate, const Vector<double>& solution, Direction& direction) const override;

      [[nodiscard]] double push_variable_to_interior(double variable_value, double lower_bound, double upper_bound) const;
      [[nodiscard]] double dual_regularization_factor() const override;
      [[nodiscard]] double compute_barrier_term_directional_derivative(const Iterate& current_iterate,
         const Vector<double>& primal_direction) const;
      void postprocess_iterate(Iterate& iterate) const override;

      // progress measures
      void set_auxiliary_measure(Iterate& iterate) const override;
      [[nodiscard]] double compute_predicted_auxiliary_reduction(const Iterate& current_iterate,
         const Vector<double>& primal_direction, double step_length) const override;

   protected:
      const OptimizationProblem& inner;
      const Parameterization& parameterization;
      const InteriorPointParameters& parameters;
      const Vector<size_t> fixed_variables{};
      const ForwardRange equality_constraints;
      const ForwardRange inequality_constraints{0};

      std::vector<double> barrier_variables_lower_bounds;
      std::vector<double> barrier_variables_upper_bounds;

      void compute_bound_dual_direction(const Iterate& current_iterate, Direction& direction) const;
      [[nodiscard]] double primal_fraction_to_boundary(const Vector<double>& current_primals, const Vector<double>& primal_direction,
         double tau) const;
      [[nodiscard]] double dual_fraction_to_boundary(const Multipliers& current_multipliers, const Multipliers& direction_multipliers,
         double tau) const;
   };
} // namespace

#endif // UNO_PRIMALDUALINTERIORPOINTPROBLEM_H