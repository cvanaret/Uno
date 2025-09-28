// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_OPTIMIZATIONPROBLEM_H
#define UNO_OPTIMIZATIONPROBLEM_H

#include "linear_algebra/MatrixOrder.hpp"
#include "linear_algebra/Norm.hpp"
#include "model/Model.hpp"
#include "optimization/SolutionStatus.hpp"
#include "optimization/LagrangianGradient.hpp"

namespace uno {
   // forward declarations
   template <typename ElementType>
   class Collection;
   template <typename IndexType>
   class COOMatrix;
   class Direction;
   class HessianModel;
   class InequalityHandlingMethod;
   class Iterate;
   class Multipliers;
   class Statistics;

   class OptimizationProblem {
   public:
      explicit OptimizationProblem(const Model& model);
      OptimizationProblem(const Model& model, size_t number_variables, size_t number_constraints);
      virtual ~OptimizationProblem() = default;

      const Model& model;
      const size_t number_variables; /*!< Number of variables */
      const size_t number_constraints; /*!< Number of constraints */

      [[nodiscard]] virtual double get_objective_multiplier() const;

      // constraint evaluations
      virtual void evaluate_constraints(Iterate& iterate, Vector<double>& constraints) const;

      // dense objective gradient
      virtual void evaluate_objective_gradient(Iterate& iterate, double* objective_gradient) const;

      // sparsity patterns of Jacobian and Hessian
      [[nodiscard]] virtual size_t number_jacobian_nonzeros() const;
      [[nodiscard]] virtual bool has_curvature(const HessianModel& hessian_model) const;
      [[nodiscard]] virtual size_t number_hessian_nonzeros(const HessianModel& hessian_model) const;
      virtual void compute_constraint_jacobian_sparsity(int* row_indices, int* column_indices, int solver_indexing,
         MatrixOrder matrix_order) const;
      virtual void compute_hessian_sparsity(const HessianModel& hessian_model, int* row_indices,
         int* column_indices, int solver_indexing) const;

      // numerical evaluations of Jacobian and Hessian
      virtual void evaluate_constraint_jacobian(Iterate& iterate, double* jacobian_values) const;
      virtual void evaluate_lagrangian_gradient(LagrangianGradient<double>& lagrangian_gradient,
         const InequalityHandlingMethod& inequality_handling_method, Iterate& iterate) const;
      virtual void evaluate_lagrangian_hessian(Statistics& statistics, HessianModel& hessian_model, const Vector<double>& primal_variables,
         const Multipliers& multipliers, double* hessian_values) const;
      virtual void compute_hessian_vector_product(HessianModel& hessian_model, const double* x, const double* vector,
         const Multipliers& multipliers, double* result) const;

      [[nodiscard]] size_t get_number_original_variables() const;
      [[nodiscard]] virtual double variable_lower_bound(size_t variable_index) const;
      [[nodiscard]] virtual double variable_upper_bound(size_t variable_index) const;
      [[nodiscard]] virtual const Vector<size_t>& get_fixed_variables() const;
      [[nodiscard]] virtual const Collection<size_t>& get_primal_regularization_variables() const;

      [[nodiscard]] virtual double constraint_lower_bound(size_t constraint_index) const;
      [[nodiscard]] virtual double constraint_upper_bound(size_t constraint_index) const;
      [[nodiscard]] virtual const Collection<size_t>& get_equality_constraints() const;
      [[nodiscard]] virtual const Collection<size_t>& get_inequality_constraints() const;
      [[nodiscard]] virtual const Collection<size_t>& get_dual_regularization_constraints() const;

      virtual void assemble_primal_dual_direction(const Iterate& current_iterate, const Vector<double>& solution, Direction& direction) const;
      [[nodiscard]] virtual double dual_regularization_factor() const;

      [[nodiscard]] static double stationarity_error(const LagrangianGradient<double>& lagrangian_gradient, double objective_multiplier,
         Norm residual_norm);
      [[nodiscard]] virtual double complementarity_error(const Vector<double>& primals, const Vector<double>& constraints,
         const Multipliers& multipliers, double shift_value, Norm residual_norm) const;

      [[nodiscard]] SolutionStatus check_first_order_convergence(const Iterate& current_iterate, double primal_tolerance,
         double dual_tolerance) const;

   protected:
      const ForwardRange primal_regularization_variables;
      const ForwardRange dual_regularization_constraints;
   };
} // namespace

#endif // UNO_OPTIMIZATIONPROBLEM_H