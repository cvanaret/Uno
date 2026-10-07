// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "NoInertiaCorrection.hpp"
#include "ingredients/subproblem/Subproblem.hpp"

namespace uno {
   void NoInertiaCorrection::initialize_statistics(Statistics& /*statistics*/) {
      // do nothing
   }

   void NoInertiaCorrection::regularize_hessian(Statistics& /*statistics*/, const Subproblem& /*subproblem*/,
         const Inertia& /*expected_inertia*/, View<double> /*hessian_values*/) {
      // do nothing
   }

   void NoInertiaCorrection::regularize_hessian(Statistics& /*statistics*/, const Subproblem& /*subproblem*/,
         const Inertia& /*expected_inertia*/, DirectSymmetricIndefiniteLinearSolver<double>& /*linear_solver*/,
         View<double> /*hessian_values*/) {
      // do nothing
   }

   void NoInertiaCorrection::regularize_augmented_matrix(Statistics& /*statistics*/, const Subproblem& /*subproblem*/,
         double /*dual_regularization_parameter*/, const Inertia& /*expected_inertia*/, BlockAugmentedMatrix block_augmented_matrix,
         Range<> /*primal_regularization_variables*/) {
      block_augmented_matrix.primal_inertia_correction.fill(0.);
      block_augmented_matrix.dual_inertia_correction.fill(0.);
   }

   void NoInertiaCorrection::regularize_augmented_matrix(Statistics& /*statistics*/, const Subproblem& /*subproblem*/,
         double /*dual_regularization_parameter*/, const Inertia& /*expected_inertia*/, DirectSymmetricIndefiniteLinearSolver<double>& /*linear_solver*/,
         BlockAugmentedMatrix block_augmented_matrix, Range<> /*primal_regularization_variables*/) {
      block_augmented_matrix.primal_inertia_correction.fill(0.);
      block_augmented_matrix.dual_inertia_correction.fill(0.);
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