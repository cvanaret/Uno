// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRIMALINERTIACORRECTION_H
#define UNO_PRIMALINERTIACORRECTION_H

#include <memory>
#include <string>
#include "InertiaCorrectionStrategy.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"

namespace uno {
   // forward declaration
   template <typename ElementType>
   class DirectSymmetricIndefiniteLinearSolver;

   class PrimalInertiaCorrection: public InertiaCorrectionStrategy {
   public:
      explicit PrimalInertiaCorrection(const Options& options);

      void initialize_statistics(Statistics& statistics) override;

      void regularize_hessian(Statistics& statistics, const Subproblem& subproblem, const Inertia& expected_inertia,
         double* hessian_values) override;
      void regularize_hessian(Statistics& statistics, const Subproblem& subproblem, const Inertia& expected_inertia,
         DirectSymmetricIndefiniteLinearSolver<double>& linear_solver, double* hessian_values) override;
      void regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         double dual_regularization_parameter, const Inertia& expected_inertia, View<double> primal_inertia_correction_block,
         View<double> dual_inertia_correction_block) override;
      void regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         double dual_regularization_parameter, const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         View<double> primal_inertia_correction_block, View<double> dual_inertia_correction_block) override;

      [[nodiscard]] bool performs_primal_regularization() const override;
      [[nodiscard]] bool performs_dual_regularization() const override;
      [[nodiscard]] double get_primal_regularization_factor() const override;
      [[nodiscard]] std::string get_name() const override;

   protected:
      const std::string& optional_linear_solver_name;
      const std::string libhsl_path;
      std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<double>> optional_linear_solver{};
      double regularization_factor{0.};
      const double regularization_initial_factor{};
      const double regularization_increase_factor{};
      const double regularization_failure_threshold{};
   };
} // namespace

#endif // UNO_PRIMALINERTIACORRECTION_H