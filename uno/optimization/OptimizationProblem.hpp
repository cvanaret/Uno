// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_OPTIMIZATIONPROBLEM_H
#define UNO_OPTIMIZATIONPROBLEM_H

#include <vector>
#include "linear_algebra/Norm.hpp"
#include "model/Model.hpp"
#include "optimization/IterateStatus.hpp"
#include "optimization/LagrangianGradient.hpp"

namespace uno {
   // forward declarations
   template <typename ElementType>
   class Collection;
   class Direction;
   class HessianModel;
   class Iterate;
   class Multipliers;
   template <typename ElementType>
   class RectangularMatrix;
   template <typename ElementType>
   class SparseVector;
   class Statistics;
   template <typename IndexType, typename ElementType>
   class SymmetricMatrix;

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
      virtual void evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const;

      // dense objective gradient
      virtual void evaluate_objective_gradient(Iterate& iterate, Vector<double>& objective_gradient) const;

      // structures of Jacobian and Hessian
      virtual void compute_hessian_structure(const HessianModel& hessian_model, Vector<size_t>& row_indices,
         Vector<size_t>& column_indices) const;

      // numerical evaluations of Jacobian and Hessian
      virtual void evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const;
      virtual void evaluate_lagrangian_gradient(LagrangianGradient<double>& lagrangian_gradient, Iterate& iterate,
         const Multipliers& multipliers) const;
      virtual void evaluate_lagrangian_hessian(Statistics& statistics, HessianModel& hessian_model, const Vector<double>& primal_variables,
         const Multipliers& multipliers, SymmetricMatrix<size_t, double>& hessian) const;
      virtual void compute_hessian_vector_product(HessianModel& hessian_model, const double* vector,
         const Multipliers& multipliers, double* result) const;

      [[nodiscard]] size_t get_number_original_variables() const;
      [[nodiscard]] virtual double variable_lower_bound(size_t variable_index) const;
      [[nodiscard]] virtual double variable_upper_bound(size_t variable_index) const;
      [[nodiscard]] virtual const Collection<size_t>& get_lower_bounded_variables() const;
      [[nodiscard]] virtual const Collection<size_t>& get_upper_bounded_variables() const;
      [[nodiscard]] virtual const Collection<size_t>& get_single_lower_bounded_variables() const;
      [[nodiscard]] virtual const Collection<size_t>& get_single_upper_bounded_variables() const;
      [[nodiscard]] virtual const Vector<size_t>& get_fixed_variables() const;
      [[nodiscard]] virtual const Collection<size_t>& get_primal_regularization_variables() const;

      [[nodiscard]] virtual double constraint_lower_bound(size_t constraint_index) const;
      [[nodiscard]] virtual double constraint_upper_bound(size_t constraint_index) const;
      [[nodiscard]] virtual const Collection<size_t>& get_equality_constraints() const;
      [[nodiscard]] virtual const Collection<size_t>& get_inequality_constraints() const;
      [[nodiscard]] virtual const Collection<size_t>& get_dual_regularization_constraints() const;

      [[nodiscard]] virtual size_t number_jacobian_nonzeros() const;
      [[nodiscard]] virtual bool has_curvature(const HessianModel& hessian_model) const;
      [[nodiscard]] virtual size_t number_hessian_nonzeros(const HessianModel& hessian_model) const;

      virtual void assemble_primal_dual_direction(const Iterate& current_iterate, const Vector<double>& solution, Direction& direction) const;
      [[nodiscard]] virtual double dual_regularization_factor() const;

      [[nodiscard]] static double stationarity_error(const LagrangianGradient<double>& lagrangian_gradient, double objective_multiplier,
         Norm residual_norm);
      [[nodiscard]] virtual double complementarity_error(const Vector<double>& primals, const std::vector<double>& constraints,
         const Multipliers& multipliers, double shift_value, Norm residual_norm) const;

      [[nodiscard]] IterateStatus check_first_order_convergence(const Iterate& current_iterate, double primal_tolerance,
         double dual_tolerance) const;

   protected:
      const ForwardRange primal_regularization_variables;
      const ForwardRange dual_regularization_constraints;
   };
} // namespace

#endif // UNO_OPTIMIZATIONPROBLEM_H