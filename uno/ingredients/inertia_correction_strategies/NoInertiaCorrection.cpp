// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "NoInertiaCorrection.hpp"
#include "ingredients/subproblem/Subproblem.hpp"

namespace uno {
   void NoInertiaCorrection::initialize_statistics(Statistics& /*statistics*/) {
      // do nothing
   }

   void NoInertiaCorrection::regularize_hessian(Statistics& /*statistics*/, const Subproblem& /*subproblem*/, const double* /*hessian_values*/,
         const Inertia& /*expected_inertia*/, double* /*primal_regularization_values*/) {
      // do nothing
   }

   void NoInertiaCorrection::regularize_hessian(Statistics& /*statistics*/, const Subproblem& /*subproblem*/, const double* /*hessian_values*/,
         const Inertia& /*expected_inertia*/, DirectSymmetricIndefiniteLinearSolver<double>& /*linear_solver*/,
         double* /*primal_regularization_values*/) {
      // do nothing
   }

   void NoInertiaCorrection::regularize_augmented_matrix(Statistics& /*statistics*/, const Subproblem& /*subproblem*/,
         const double* /*augmented_matrix_values*/, double /*dual_regularization_parameter*/,
         const Inertia& /*expected_inertia*/, double* /*primal_regularization_values*/,
         double* /*dual_regularization_values*/) {
      // do nothing
   }

   void NoInertiaCorrection::regularize_augmented_matrix(Statistics& /*statistics*/, const Subproblem& subproblem,
         const double* /*augmented_matrix_values*/, double /*dual_regularization_parameter*/,
         const Inertia& /*expected_inertia*/, DirectSymmetricIndefiniteLinearSolver<double>& /*linear_solver*/,
         double* primal_regularization_values, double* dual_regularization_values) {
      for (size_t index: Range(subproblem.get_primal_regularization_variables().size())) {
         primal_regularization_values[index] = 0.;
      }
      for (size_t index: Range(subproblem.get_dual_regularization_constraints().size())) {
         dual_regularization_values[index] = 0.;
      }
   }

   [[nodiscard]] bool NoInertiaCorrection::performs_primal_regularization() const {
      return false;
   }

   [[nodiscard]] bool NoInertiaCorrection::performs_dual_regularization() const {
      return false;
   }

   [[nodiscard]] double NoInertiaCorrection::get_primal_regularization_factor() const {
      return 0.;
   }

   [[nodiscard]] std::string NoInertiaCorrection::get_name() const {
      return "no";
   }
} // namespace