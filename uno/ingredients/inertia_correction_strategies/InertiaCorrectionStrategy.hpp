// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INERTIACORRECTIONSTRATEGY_H
#define UNO_INERTIACORRECTIONSTRATEGY_H

#include <string>
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "linear_algebra/View.hpp"
#include "symbolic/Range.hpp"

namespace uno {
   // forward declarations
   class Inertia;
   class Options;
   class Statistics;

   class InertiaCorrectionStrategy {
   public:
      InertiaCorrectionStrategy() = default;
      virtual ~InertiaCorrectionStrategy() = default;

      virtual void initialize_statistics(Statistics& statistics) = 0;

      virtual void regularize_hessian(Statistics& statistics, const Subproblem& subproblem, const Inertia& expected_inertia,
         View<double> hessian_values) = 0;
      virtual void regularize_hessian(Statistics& statistics, const Subproblem& subproblem, const Inertia& expected_inertia,
         DirectSymmetricIndefiniteLinearSolver<double>& linear_solver, View<double> primal_regularization_values) = 0;
      virtual void regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         double dual_regularization_parameter, const Inertia& expected_inertia, BlockAugmentedMatrix block_augmented_matrix,
         Range<> primal_regularization_variables) = 0;
      virtual void regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         double dual_regularization_parameter, const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
        BlockAugmentedMatrix block_augmented_matrix, Range<> primal_regularization_variables) = 0;

      // the augmented matrix factorized with the current regularization is deemed singular by the caller (e.g. iterative
      // refinement failed): increase the regularization and refactorize the matrix until its inertia is correct.
      // Returns false if the strategy cannot handle singular matrices
      [[nodiscard]] virtual bool regularize_singular_augmented_matrix(Statistics& /*statistics*/, const Subproblem& /*subproblem*/,
            double /*dual_regularization_parameter*/, const Inertia& /*expected_inertia*/,
            DirectSymmetricIndefiniteLinearSolver<double>& /*linear_solver*/, BlockAugmentedMatrix /*block_augmented_matrix*/,
            Range<> /*primal_regularization_variables*/) {
         return false;
      }

      // the augmented matrix must be refactorized with the current regularization factors (e.g. after the pivoting
      // tolerance of the linear solver was increased). Strategies that correct the inertia check it again and increase
      // the regularization if needed; by default, the matrix is simply refactorized
      virtual void refactorize_augmented_matrix(Statistics& /*statistics*/, const Subproblem& /*subproblem*/,
            double /*dual_regularization_parameter*/, const Inertia& /*expected_inertia*/,
            DirectSymmetricIndefiniteLinearSolver<double>& linear_solver, BlockAugmentedMatrix /*block_augmented_matrix*/,
            Range<> /*primal_regularization_variables*/) {
         linear_solver.do_numerical_factorization(false);
      }

      [[nodiscard]] virtual bool performs_primal_regularization() const = 0;
      [[nodiscard]] virtual bool performs_dual_regularization() const = 0;
      [[nodiscard]] virtual double get_primal_regularization_factor() const = 0;
      [[nodiscard]] virtual std::string get_name() const = 0;
   };
} // namespace

#endif // UNO_INERTIACORRECTIONSTRATEGY_H