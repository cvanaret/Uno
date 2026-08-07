// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INERTIACORRECTIONSTRATEGY_H
#define UNO_INERTIACORRECTIONSTRATEGY_H

#include <string>
#include "linear_algebra/View.hpp"

namespace uno {
   // forward declarations
   template <typename ElementType>
   class DirectSymmetricIndefiniteLinearSolver;
   class Inertia;
   class Options;
   class Statistics;
   class Subproblem;

   class InertiaCorrectionStrategy {
   public:
      InertiaCorrectionStrategy() = default;
      virtual ~InertiaCorrectionStrategy() = default;

      virtual void initialize_statistics(Statistics& statistics) = 0;

      virtual void regularize_hessian(Statistics& statistics, const Subproblem& subproblem, const Inertia& expected_inertia,
         View<double> hessian_values) = 0;
      virtual void regularize_hessian(Statistics& statistics, const Subproblem& subproblem, const Inertia& expected_inertia,
         DirectSymmetricIndefiniteLinearSolver<double>& linear_solver, View<double> hessian_values) = 0;
      virtual void regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         double dual_regularization_parameter, const Inertia& expected_inertia, View<double> primal_inertia_correction_block,
         View<double> dual_inertia_correction_block) = 0;
      virtual void regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         double dual_regularization_parameter, const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
        View<double> primal_inertia_correction_block, View<double> dual_inertia_correction_block) = 0;

      [[nodiscard]] virtual bool performs_primal_regularization() const = 0;
      [[nodiscard]] virtual bool performs_dual_regularization() const = 0;
      [[nodiscard]] virtual double get_primal_regularization_factor() const = 0;
      [[nodiscard]] virtual std::string get_name() const = 0;
   };
} // namespace

#endif // UNO_INERTIACORRECTIONSTRATEGY_H