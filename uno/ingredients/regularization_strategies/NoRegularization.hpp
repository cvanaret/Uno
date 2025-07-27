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

      void initialize_statistics(Statistics& /*statistics*/, const Options& /*options*/) override {
         // do nothing
      }

      void regularize_hessian(Statistics& /*statistics*/, const Subproblem& /*subproblem*/, const Vector<double>& /*hessian_values*/,
            const Inertia& /*expected_inertia*/) override {
         // do nothing
      }

      void regularize_hessian(Statistics& /*statistics*/, const Subproblem& /*subproblem*/, const Vector<double>& /*hessian_values*/,
            const Inertia& /*expected_inertia*/, DirectSymmetricIndefiniteLinearSolver<size_t, double>& /*linear_solver*/) override {
         // do nothing
      }

      void regularize_augmented_matrix(Statistics& /*statistics*/, const Subproblem& /*subproblem*/,
            const Vector<double>& /*augmented_matrix_values*/, ElementType /*dual_regularization_parameter*/,
            const Inertia& /*expected_inertia*/, VectorView<Vector<double>&> /*primal_regularization_values*/,
            VectorView<Vector<double>&> /*dual_regularization_values*/) override {
         // do nothing
      }

      void regularize_augmented_matrix(Statistics& /*statistics*/, const Subproblem& subproblem,
            const Vector<double>& /*augmented_matrix_values*/, ElementType /*dual_regularization_parameter*/,
            const Inertia& /*expected_inertia*/, DirectSymmetricIndefiniteLinearSolver<size_t, double>& /*linear_solver*/,
            VectorView<Vector<double>&> primal_regularization_values, VectorView<Vector<double>&> dual_regularization_values) override {
         for (size_t index: Range(subproblem.get_primal_regularization_variables().size())) {
            primal_regularization_values[index] = 0.;
         }
         for (size_t index: Range(subproblem.get_dual_regularization_constraints().size())) {
            dual_regularization_values[index] = 0.;
         }
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