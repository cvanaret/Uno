// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRIMALDUALINERTIACORRECTION_H
#define UNO_PRIMALDUALINERTIACORRECTION_H

#include <string>
#include "InertiaCorrectionStrategy.hpp"
#include "UnstableRegularization.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"

namespace uno {
   class PrimalDualInertiaCorrection: public InertiaCorrectionStrategy {
   public:
      explicit PrimalDualInertiaCorrection(const Options& options);

      void initialize_statistics(Statistics& statistics) override;

      void regularize_hessian(Statistics& statistics, const Subproblem& subproblem, const double* hessian_values,
         const Inertia& expected_inertia, double* primal_regularization_values) override;
      void regularize_hessian(Statistics& statistics, const Subproblem& subproblem, const double* hessian_values,
         const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         double* primal_regularization_values) override;
      void regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         const double* augmented_matrix_values, double dual_regularization_parameter,
         const Inertia& expected_inertia, double* primal_regularization_values,
         double* dual_regularization_values) override;
      void regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         const double* augmented_matrix_values, double dual_regularization_parameter,
         const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         double* primal_regularization_values, double* dual_regularization_values) override;

      [[nodiscard]] bool performs_primal_regularization() const override;
      [[nodiscard]] bool performs_dual_regularization() const override;
      [[nodiscard]] double get_primal_regularization_factor() const override;
      [[nodiscard]] std::string get_name() const override;

   protected:
      const std::string& optional_linear_solver_name;
      std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<double>> optional_linear_solver{};
      double primal_regularization{0.};
      double dual_regularization{0.};
      double previous_primal_regularization{0.};
      const double regularization_failure_threshold;
      const double primal_regularization_initial_factor;
      const double dual_regularization_fraction;
      const double primal_regularization_lb;
      const double primal_regularization_decrease_factor;
      const double primal_regularization_fast_increase_factor;
      const double primal_regularization_slow_increase_factor;
      const size_t threshold_unsuccessful_attempts;
   };
} // namespace

#endif // UNO_PRIMALDUALINERTIACORRECTION_H