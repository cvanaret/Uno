// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_OPTIMIZATIONPROBLEM_H
#define UNO_OPTIMIZATIONPROBLEM_H

#include <memory>
#include <vector>
#include "ingredients/inertia_correction_strategies/Inertia.hpp"
#include "linear_algebra/MatrixOrder.hpp"
#include "linear_algebra/Norm.hpp"
#include "optimization/SolutionStatus.hpp"
#include "../interfaces/C/uno_int.h"
#include "model/Model.hpp"

namespace uno {
   // forward declarations
   template <typename ElementType>
   class Collection;
   class Direction;
   class Evaluations;
   class SolverWorkspace;
   class HessianModel;
   class Iterate;
   class Model;
   class Multipliers;
   class Parameterization;
   class Statistics;
   template <typename ElementType>
   class Vector;

   class OptimizationProblem {
   public:
      explicit OptimizationProblem(const Model& model);
      OptimizationProblem(const Model& model, size_t number_variables, size_t number_constraints);
      virtual ~OptimizationProblem() = default;
      virtual std::unique_ptr<OptimizationProblem> clone() const;

      const Model& model;
      const size_t number_variables; /*!< Number of variables */
      const size_t number_constraints; /*!< Number of constraints */

      [[nodiscard]] virtual double get_objective_multiplier() const;
      [[nodiscard]] virtual bool has_inequality_constraints() const;
      [[nodiscard]] virtual bool has_bound_constraints() const;

      virtual void generate_initial_iterate(Iterate& initial_iterate, Evaluations& evaluations) const;
      virtual void postprocess_iterate(Iterate& iterate) const;

      // sparsity patterns of Jacobian and Hessian
      [[nodiscard]] virtual size_t number_jacobian_nonzeros() const;
      [[nodiscard]] virtual bool has_curvature(const HessianModel& hessian_model) const;
      [[nodiscard]] virtual size_t number_hessian_nonzeros(const HessianModel& hessian_model) const;
      virtual void compute_jacobian_sparsity(uno_int* row_indices, uno_int* column_indices, uno_int row_offset,
         uno_int column_offset, uno_int solver_indexing, MatrixOrder matrix_order) const;
      virtual void compute_hessian_sparsity(const HessianModel& hessian_model, uno_int* row_indices,
         uno_int* column_indices, uno_int solver_indexing) const;

      // numerical evaluations of constraints, objective gradient, Jacobian and Hessian
      virtual void evaluate_constraints(const Iterate& iterate, double* constraints, Evaluations& evaluations) const;
      virtual void evaluate_objective_gradient(const Iterate& iterate, double* objective_gradient, Evaluations& evaluations) const;
      virtual void evaluate_jacobian(const Vector<double>& primals, double* jacobian_values, Evaluations& evaluations) const;
      virtual void evaluate_lagrangian_gradient(const Iterate& iterate, Evaluations& evaluations,
         Vector<double>& lagrangian_gradient) const;
      virtual void evaluate_lagrangian_hessian(Statistics& statistics, HessianModel& hessian_model,
         const Vector<double>& primal_variables, const Multipliers& multipliers, double* hessian_values) const;

      // linear operators
      virtual void compute_jacobian_vector_product(const double* vector, double* result, const Evaluations& evaluations) const;
      virtual void compute_jacobian_transposed_vector_product(const double* vector, double* result,
         const Evaluations& evaluations) const;
      virtual void compute_hessian_vector_product(HessianModel& hessian_model, const double* x, const double* vector,
         const Multipliers& multipliers, double* result) const;

      [[nodiscard]] size_t get_number_original_variables() const;
      [[nodiscard]] virtual const std::vector<double>& get_variables_lower_bounds() const;
      [[nodiscard]] virtual const std::vector<double>& get_variables_upper_bounds() const;
      [[nodiscard]] virtual const Vector<size_t>& get_fixed_variables() const;
      [[nodiscard]] virtual const Collection<size_t>& get_primal_regularization_variables() const;

      [[nodiscard]] virtual const std::vector<double>& get_constraints_lower_bounds() const;
      [[nodiscard]] virtual const std::vector<double>& get_constraints_upper_bounds() const;
      [[nodiscard]] virtual const Collection<size_t>& get_equality_constraints() const;
      [[nodiscard]] virtual const Collection<size_t>& get_inequality_constraints() const;
      [[nodiscard]] virtual const Collection<size_t>& get_dual_regularization_constraints() const;

      [[nodiscard]] virtual Inertia get_inertia() const;

      virtual void assemble_primal_dual_direction(const Iterate& current_iterate, const Vector<double>& solution, Direction& direction) const;
      [[nodiscard]] virtual double dual_regularization_factor() const;

      [[nodiscard]] virtual double complementarity_error(const Vector<double>& primals, const Vector<double>& constraints,
         const Multipliers& multipliers, double shift_value, Norm residual_norm) const;
      [[nodiscard]] virtual double compute_centrality_error(const Vector<double>& primals, const Multipliers& multipliers,
         double shift) const;

      [[nodiscard]] virtual SolutionStatus check_first_order_convergence(const Iterate& current_iterate, double primal_tolerance,
         double dual_tolerance) const;

      // progress measures
      virtual void set_infeasibility_measure(Iterate& iterate, Evaluations& evaluations, Norm norm) const;
      virtual void set_objective_measure(Iterate& iterate, Evaluations& evaluations) const;
      virtual void set_auxiliary_measure(Iterate& iterate) const;
      [[nodiscard]] virtual double compute_predicted_auxiliary_reduction(const Iterate& current_iterate,
         const Vector<double>& primal_direction, double step_length) const;

   protected:
      const ForwardRange primal_regularization_variables;
      const ForwardRange dual_regularization_constraints;
   };
} // namespace

#endif // UNO_OPTIMIZATIONPROBLEM_H