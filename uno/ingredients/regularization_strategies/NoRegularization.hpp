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
      explicit NoRegularization() = default;

      void initialize_statistics(Statistics& statistics, const Options& options) override;

      void regularize_hessian(Statistics& statistics, DirectSymmetricIndefiniteLinearSolver<size_t, ElementType>& linear_solver,
         SymmetricMatrix<size_t, ElementType>& hessian) override;
      void regularize_augmented_matrix(Statistics& statistics, DirectSymmetricIndefiniteLinearSolver<size_t, ElementType>& linear_solver,
         SymmetricMatrix<size_t, ElementType>& augmented_matrix, size_t size_primal_block, size_t size_dual_block,
         ElementType dual_regularization_parameter) override;
   };

   template <typename ElementType>
   void NoRegularization<ElementType>::initialize_statistics(Statistics& /*statistics*/, const Options& /*options*/) { }

   template <typename ElementType>
   void NoRegularization<ElementType>::regularize_hessian(Statistics& /*statistics*/,
         DirectSymmetricIndefiniteLinearSolver<size_t, ElementType>& /*linear_solver*/, SymmetricMatrix<size_t, ElementType>& /*hessian*/) {
      // do nothing
   }

   template <typename ElementType>
   void NoRegularization<ElementType>::regularize_augmented_matrix(Statistics& /*statistics*/,
         DirectSymmetricIndefiniteLinearSolver<size_t, ElementType>& /*linear_solver*/, SymmetricMatrix<size_t, ElementType>& /*augmented_matrix*/,
         size_t /*size_primal_block*/, size_t /*size_dual_block*/, ElementType /*dual_regularization_parameter*/) {
      // do nothing
   }
} // namespace

#endif // UNO_NOREGULARIZATION_H
