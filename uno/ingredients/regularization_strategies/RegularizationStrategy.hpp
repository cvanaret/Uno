// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_REGULARIZATIONSTRATEGY_H
#define UNO_REGULARIZATIONSTRATEGY_H

namespace uno {
   // forward declarations
   template <typename IndexType, typename ElementType>
   class DirectSymmetricIndefiniteLinearSolver;
   class Statistics;
   template <typename IndexType, typename ElementType>
   class SymmetricMatrix;

   template <typename ElementType>
   class RegularizationStrategy {
   public:
      RegularizationStrategy() = default;
      virtual ~RegularizationStrategy();

      virtual void regularize_matrix(Statistics& statistics, DirectSymmetricIndefiniteLinearSolver<size_t, ElementType>& linear_solver,
            SymmetricMatrix<size_t, ElementType>& matrix, size_t size_primal_block, size_t size_dual_block,
            ElementType dual_regularization_parameter) = 0;
   };

   template <typename ElementType>
   RegularizationStrategy<ElementType>::~RegularizationStrategy<ElementType>() { }
} // namespace

#endif // UNO_REGULARIZATIONSTRATEGY_H