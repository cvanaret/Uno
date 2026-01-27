// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_NOINERTIACORRECTION_H
#define UNO_NOINERTIACORRECTION_H

#include "InertiaCorrectionStrategy.hpp"

namespace uno {
   class NoInertiaCorrection: public InertiaCorrectionStrategy {
   public:
      explicit NoInertiaCorrection() = default;

      void initialize_statistics(Statistics& /*statistics*/) override;

      void regularize_hessian(Statistics& /*statistics*/, const Subproblem& /*subproblem*/, const double* /*hessian_values*/,
            const Inertia& /*expected_inertia*/, double* /*primal_regularization_values*/) override;
      void regularize_hessian(Statistics& /*statistics*/, const Subproblem& /*subproblem*/, const double* /*hessian_values*/,
            const Inertia& /*expected_inertia*/, DirectSymmetricIndefiniteLinearSolver<double>& /*linear_solver*/,
            double* /*primal_regularization_values*/) override;
      void regularize_augmented_matrix(Statistics& /*statistics*/, const Subproblem& /*subproblem*/,
            const double* /*augmented_matrix_values*/, double /*dual_regularization_parameter*/,
            const Inertia& /*expected_inertia*/, double* /*primal_regularization_values*/,
            double* /*dual_regularization_values*/) override;
      void regularize_augmented_matrix(Statistics& /*statistics*/, const Subproblem& subproblem,
            const double* /*augmented_matrix_values*/, double /*dual_regularization_parameter*/,
            const Inertia& /*expected_inertia*/, DirectSymmetricIndefiniteLinearSolver<double>& /*linear_solver*/,
            double* primal_regularization_values, double* dual_regularization_values) override;

      [[nodiscard]] bool performs_primal_regularization() const override;
      [[nodiscard]] bool performs_dual_regularization() const override;
      [[nodiscard]] double get_primal_regularization_factor() const override;
      [[nodiscard]] std::string get_name() const override;
   };
} // namespace

#endif // UNO_NOINERTIACORRECTION_H