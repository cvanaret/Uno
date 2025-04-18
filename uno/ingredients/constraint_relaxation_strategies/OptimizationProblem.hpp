// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_OPTIMIZATIONPROBLEM_H
#define UNO_OPTIMIZATIONPROBLEM_H

#include <vector>
#include "linear_algebra/Norm.hpp"
#include "model/Model.hpp"
#include "optimization/LagrangianGradient.hpp"
#include "symbolic/Expression.hpp"

namespace uno {
   // forward declarations
   template <typename ElementType>
   class Collection;
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
      OptimizationProblem(const Model& model, size_t number_variables, size_t number_constraints);
      virtual ~OptimizationProblem() = default;

      const Model& model;
      const size_t number_variables; /*!< Number of variables */
      const size_t number_constraints; /*!< Number of constraints */

      [[nodiscard]] bool is_constrained() const;
      [[nodiscard]] bool has_inequality_constraints() const;
      [[nodiscard]] bool has_fixed_variables() const;

      // function evaluations
      [[nodiscard]] virtual double get_objective_multiplier() const;
      virtual void evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const;
      virtual void evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const;
      virtual void evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const;
      virtual void evaluate_lagrangian_hessian(Statistics& statistics, HessianModel& hessian_model, const Vector<double>& primal_variables,
         const Multipliers& multipliers, SymmetricMatrix<size_t, double>& hessian) const;
      virtual void compute_hessian_vector_product(HessianModel& hessian_model, const Vector<double>& vector,
         const Multipliers& multipliers, Vector<double>& result) const;

      [[nodiscard]] size_t get_number_original_variables() const;
      [[nodiscard]] virtual double variable_lower_bound(size_t variable_index) const;
      [[nodiscard]] virtual double variable_upper_bound(size_t variable_index) const;
      [[nodiscard]] virtual const Collection<size_t>& get_lower_bounded_variables() const;
      [[nodiscard]] virtual const Collection<size_t>& get_upper_bounded_variables() const;
      [[nodiscard]] virtual const Collection<size_t>& get_single_lower_bounded_variables() const;
      [[nodiscard]] virtual const Collection<size_t>& get_single_upper_bounded_variables() const;

      [[nodiscard]] virtual double constraint_lower_bound(size_t constraint_index) const;
      [[nodiscard]] virtual double constraint_upper_bound(size_t constraint_index) const;
      [[nodiscard]] virtual const Collection<size_t>& get_inequality_constraints() const;

      [[nodiscard]] virtual size_t number_objective_gradient_nonzeros() const;
      [[nodiscard]] virtual size_t number_jacobian_nonzeros() const;
      [[nodiscard]] virtual size_t number_hessian_nonzeros() const;

      [[nodiscard]] static double stationarity_error(const LagrangianGradient<double>& lagrangian_gradient, double objective_multiplier,
            Norm residual_norm);
      virtual void evaluate_lagrangian_gradient(LagrangianGradient<double>& lagrangian_gradient, Iterate& iterate, const Multipliers& multipliers) const;
      [[nodiscard]] virtual double complementarity_error(const Vector<double>& primals, const std::vector<double>& constraints,
            const Multipliers& multipliers, double shift_value, Norm residual_norm) const;
   };
} // namespace

#endif // UNO_OPTIMIZATIONPROBLEM_H