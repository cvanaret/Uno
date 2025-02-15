// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_NOREGULARIZATION_H
#define UNO_NOREGULARIZATION_H

#include <cassert>
#include "RegularizationStrategy.hpp"

namespace uno {
   template <typename ElementType>
   class NoRegularization: public RegularizationStrategy<ElementType> {
   public:
      explicit NoRegularization();
      void regularize_matrix(Statistics& statistics, DirectSymmetricIndefiniteLinearSolver<size_t, ElementType>& linear_solver,
            SymmetricMatrix<size_t, ElementType>& matrix, size_t size_primal_block, size_t size_dual_block,
            ElementType dual_regularization_parameter) override;
   };

   template <typename ElementType>
   NoRegularization<ElementType>::NoRegularization():
         RegularizationStrategy<ElementType>() {
   }

   template <typename ElementType>
   void NoRegularization<ElementType>::regularize_matrix(Statistics& /*statistics*/,
         DirectSymmetricIndefiniteLinearSolver<size_t, ElementType>& /*linear_solver*/, SymmetricMatrix<size_t, ElementType>& /*matrix*/,
         size_t /*size_primal_block*/, size_t /*size_dual_block*/, ElementType /*dual_regularization_parameter*/) {
      // do nothing
   }
} // namespace

#endif // UNO_NOREGULARIZATION_H
