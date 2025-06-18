// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_REGULARIZATIONSTRATEGY_H
#define UNO_REGULARIZATIONSTRATEGY_H

#include "Inertia.hpp"

namespace uno {
   // forward declarations
   template <typename ElementType>
   class Collection;
   template <typename IndexType, typename ElementType>
   class DirectSymmetricIndefiniteLinearSolver;
   class HessianModel;
   class OptimizationProblem;
   class Options;
   class Statistics;
   template <typename IndexType, typename ElementType>
   class SymmetricMatrix;

   template <typename ElementType>
   class RegularizationStrategy {
   public:
      RegularizationStrategy() = default;
      virtual ~RegularizationStrategy() = default;

      virtual void initialize_memory(const OptimizationProblem& problem, const HessianModel& hessian_model) = 0;
      virtual void initialize_statistics(Statistics& statistics, const Options& options) = 0;

      virtual void regularize_hessian(Statistics& statistics, SymmetricMatrix<size_t, ElementType>& hessian,
         const Collection<size_t>& indices, const Inertia& expected_inertia) = 0;
      virtual void regularize_hessian(Statistics& statistics, SymmetricMatrix<size_t, ElementType>& hessian,
         const Collection<size_t>& indices, const Inertia& expected_inertia,
         DirectSymmetricIndefiniteLinearSolver<size_t, double>& linear_solver) = 0;
      virtual void regularize_augmented_matrix(Statistics& statistics, SymmetricMatrix<size_t, ElementType>& augmented_matrix,
         const Collection<size_t>& primal_indices, const Collection<size_t>& dual_indices,
         ElementType dual_regularization_parameter, const Inertia& expected_inertia) = 0;
      virtual void regularize_augmented_matrix(Statistics& statistics, SymmetricMatrix<size_t, ElementType>& augmented_matrix,
         const Collection<size_t>& primal_block, const Collection<size_t>& dual_block,
         ElementType dual_regularization_parameter, const Inertia& expected_inertia,
         DirectSymmetricIndefiniteLinearSolver<size_t, double>& linear_solver) = 0;

      [[nodiscard]] virtual bool performs_primal_regularization() const = 0;
      [[nodiscard]] virtual bool performs_dual_regularization() const = 0;
      [[nodiscard]] virtual std::string get_name() const = 0;
   };
} // namespace

#endif // UNO_REGULARIZATIONSTRATEGY_H