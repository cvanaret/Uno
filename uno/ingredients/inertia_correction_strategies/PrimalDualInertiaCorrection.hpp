// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRIMALDUALINERTIACORRECTION_H
#define UNO_PRIMALDUALINERTIACORRECTION_H

#include <string>
#include "InertiaCorrectionStrategy.hpp"

namespace uno {
   // forward declaration
   template <typename ElementType>
   class DirectSymmetricIndefiniteLinearSolver;

   class PrimalDualInertiaCorrection: public InertiaCorrectionStrategy {
   public:
      explicit PrimalDualInertiaCorrection(const Options& options);

      void initialize_statistics(Statistics& statistics) override;

      void regularize_hessian(Statistics& statistics, const Subproblem& subproblem, const Inertia& expected_inertia,
         View<double> hessian_values) override;
      void regularize_hessian(Statistics& statistics, const Subproblem& subproblem, const Inertia& expected_inertia,
         DirectSymmetricIndefiniteLinearSolver<double>& linear_solver, View<double> primal_regularization_values) override;
      void regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem, double dual_regularization_parameter,
         const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         BlockAugmentedMatrix block_augmented_matrix, Range<> primal_regularization_variables) override;

      [[nodiscard]] bool regularize_singular_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         double dual_regularization_parameter, const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         BlockAugmentedMatrix block_augmented_matrix, Range<> primal_regularization_variables) override;

      void refactorize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         double dual_regularization_parameter, const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         BlockAugmentedMatrix block_augmented_matrix, Range<> primal_regularization_variables) override;

      [[nodiscard]] bool performs_primal_regularization() const override;
      [[nodiscard]] bool performs_dual_regularization() const override;
      [[nodiscard]] double get_primal_regularization_factor() const override;
      [[nodiscard]] std::string get_name() const override;

   protected:
      // structural degeneracy of the Hessian and the Jacobian, detected as in IPOPT (PDPerturbationHandler)
      enum class Degeneracy {NOT_YET_DETERMINED, NOT_DEGENERATE, DEGENERATE};
      // regularization tested on the current matrix to determine the degeneracy
      enum class DegeneracyTest {NONE, NO_REGULARIZATION, DUAL_ONLY, PRIMAL_ONLY, PRIMAL_AND_DUAL};

      const Options& options;
      std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<double>> optional_linear_solver{};
      double primal_regularization{0.};
      double dual_regularization{0.};
      double previous_primal_regularization{0.};
      Degeneracy hessian_degeneracy{Degeneracy::NOT_YET_DETERMINED};
      Degeneracy jacobian_degeneracy{Degeneracy::NOT_YET_DETERMINED};
      DegeneracyTest degeneracy_test{DegeneracyTest::NONE};
      size_t number_degenerate_matrices{0};
      static constexpr size_t degeneracy_threshold{3}; // number of degenerate matrices before the degeneracy is assumed
      bool quality_improved_for_current_matrix{false};
      const double regularization_failure_threshold;
      const double primal_regularization_initial_factor;
      const double dual_regularization_fraction;
      const double primal_regularization_lb;
      const double primal_regularization_decrease_factor;
      const double primal_regularization_fast_increase_factor;
      const double primal_regularization_slow_increase_factor;
      const bool regularize_all_variables;

      void regularize_until_correct_inertia(Statistics& statistics, double dual_regularization_value,
         const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         BlockAugmentedMatrix& block_augmented_matrix, Range<> primal_regularization_variables, bool assume_singular);
      void regularize_for_singularity(double dual_regularization_value);
      void regularize_for_wrong_inertia(double dual_regularization_value);
      [[nodiscard]] bool increase_primal_regularization();
      void finalize_degeneracy_test();
   };
} // namespace

#endif // UNO_PRIMALDUALINERTIACORRECTION_H