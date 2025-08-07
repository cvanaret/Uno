// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_NOREGULARIZATION_H
#define UNO_NOREGULARIZATION_H

#include "RegularizationStrategy.hpp"

namespace uno {
   template <typename ElementType>
   class NoRegularization: public RegularizationStrategy<ElementType> {
   public:
      explicit NoRegularization() = default;

      void initialize_memory(const OptimizationProblem& /*problem*/, const HessianModel& /*hessian_model*/) override {
         // do nothing
      }

      void initialize_statistics(Statistics& /*statistics*/, const Options& /*options*/) override {
         // do nothing
      }

      void regularize_hessian(Statistics& /*statistics*/, const Subproblem& subproblem, SymmetricMatrix<size_t, ElementType>& /*hessian*/,
            const Inertia& /*expected_inertia*/) override {
         // do nothing
      }

      void regularize_hessian(Statistics& /*statistics*/, const Subproblem& subproblem, SymmetricMatrix<size_t, ElementType>& /*hessian*/,
            const Inertia& /*expected_inertia*/, DirectSymmetricIndefiniteLinearSolver<size_t, double>& /*linear_solver*/) override {
         // do nothing
      }

      void regularize_augmented_matrix(Statistics& /*statistics*/, const Subproblem& subproblem,
            SymmetricMatrix<size_t, ElementType>& /*augmented_matrix*/, ElementType /*dual_regularization_parameter*/,
            const Inertia& /*expected_inertia*/) override {
         // do nothing
      }

      void regularize_augmented_matrix(Statistics& /*statistics*/, const Subproblem& subproblem,
            SymmetricMatrix<size_t, ElementType>& /*augmented_matrix*/, ElementType /*dual_regularization_parameter*/,
            const Inertia& /*expected_inertia*/, DirectSymmetricIndefiniteLinearSolver<size_t, double>& /*linear_solver*/) override {
         // do nothing
      }

      [[nodiscard]] bool performs_primal_regularization() const override {
         return false;
      }

      [[nodiscard]] bool performs_dual_regularization() const override {
         return false;
      }

      [[nodiscard]] double get_primal_regularization_factor() const override {
         return 0.;
      }

      [[nodiscard]] std::string get_name() const override {
         return "no";
      }
   };
} // namespace

#endif // UNO_NOREGULARIZATION_H