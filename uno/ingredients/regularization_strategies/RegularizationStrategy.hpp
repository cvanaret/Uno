// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_REGULARIZATIONSTRATEGY_H
#define UNO_REGULARIZATIONSTRATEGY_H

#include "Inertia.hpp"

namespace uno {
   // forward declarations
   template <typename ElementType>
   class Collection;
   template <typename ElementType>
   class DirectSymmetricIndefiniteLinearSolver;
   class HessianModel;
   class OptimizationProblem;
   class Options;
   class Statistics;
   class Subproblem;
   template <typename ElementType>
   class Vector;

   template <typename ElementType>
   class RegularizationStrategy {
   public:
      RegularizationStrategy() = default;
      virtual ~RegularizationStrategy() = default;

      virtual void initialize_statistics(Statistics& statistics, const Options& options) = 0;

      virtual void regularize_hessian(Statistics& statistics, const Subproblem& subproblem, const double* hessian_values,
         const Inertia& expected_inertia, double* primal_regularization_values) = 0;
      virtual void regularize_hessian(Statistics& statistics, const Subproblem& subproblem, const double* hessian_values,
         const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         double* primal_regularization_values) = 0;
      virtual void regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         const double* augmented_matrix_values, ElementType dual_regularization_parameter,
         const Inertia& expected_inertia, double* primal_regularization_values,
         double* dual_regularization_values) = 0;
      virtual void regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         const double* augmented_matrix_values, ElementType dual_regularization_parameter,
         const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         double* primal_regularization_values, double* dual_regularization_values) = 0;

      [[nodiscard]] virtual bool performs_primal_regularization() const = 0;
      [[nodiscard]] virtual bool performs_dual_regularization() const = 0;
      [[nodiscard]] virtual double get_primal_regularization_factor() const = 0;

      [[nodiscard]] virtual std::string get_name() const = 0;
   };
} // namespace

#endif // UNO_REGULARIZATIONSTRATEGY_H