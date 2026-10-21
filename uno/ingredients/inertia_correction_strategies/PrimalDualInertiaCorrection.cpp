// Copyright (c) 2025-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "PrimalDualInertiaCorrection.hpp"
#include "UnstableInertiaCorrection.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem_solvers/LinearSystem.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   PrimalDualInertiaCorrection::PrimalDualInertiaCorrection(const Options& options):
         InertiaCorrectionStrategy(),
         options(options),
         regularization_failure_threshold(options.get_double("regularization_failure_threshold")),
         primal_regularization_initial_factor(options.get_double("primal_regularization_initial_factor")),
         dual_regularization_fraction(options.get_double("dual_regularization_fraction")),
         primal_regularization_lb(options.get_double("primal_regularization_lb")),
         primal_regularization_decrease_factor(options.get_double("primal_regularization_decrease_factor")),
         primal_regularization_fast_increase_factor(options.get_double("primal_regularization_fast_increase_factor")),
         primal_regularization_slow_increase_factor(options.get_double("primal_regularization_slow_increase_factor")),
         threshold_unsuccessful_attempts(options.get_unsigned_int("threshold_unsuccessful_attempts")),
         regularize_all_variables(options.get_bool("regularize_all_variables")) {
   }

   void PrimalDualInertiaCorrection::initialize_statistics(Statistics& statistics) {
      statistics.add_column("Prim reg", Statistics::double_width + 1, 2, /* is_extended = */ true);
      statistics.add_column("Dual reg", Statistics::double_width + 1, 2, /* is_extended = */ true);
   }

   // regularization of the Hessian only (e.g. for inequality-constrained QP subproblems): there is no dual block, only
   // the primal regularization is computed
   void PrimalDualInertiaCorrection::regularize_hessian(Statistics& statistics, const Subproblem& subproblem,
         const Inertia& expected_inertia, View<double> hessian_values) {
      // pick the member linear solver (system with the regularized Hessian only)
      if (this->optional_linear_solver == nullptr) {
         this->optional_linear_solver = SymmetricIndefiniteLinearSolverFactory::create(this->options);
         this->optional_linear_solver->get_linear_system().initialize_hessian(subproblem);
         this->optional_linear_solver->initialize_memory();
         this->optional_linear_solver->do_symbolic_analysis();
      }
      // copy the Hessian into the linear system
      double* matrix = this->optional_linear_solver->get_linear_system().matrix_values.data();
      const size_t number_hessian_nonzeros = subproblem.number_hessian_nonzeros();
      for (size_t nonzero_index: Range(number_hessian_nonzeros)) {
         matrix[nonzero_index] = hessian_values[nonzero_index];
      }
      // the regularization terms follow the Hessian terms in the linear system
      const size_t number_regularization_nonzeros = subproblem.number_primal_inertia_correction_nonzeros();
      View<double> primal_regularization_values(matrix + number_hessian_nonzeros, number_regularization_nonzeros);
      this->regularize_hessian(statistics, subproblem, expected_inertia, *this->optional_linear_solver, primal_regularization_values);
      // copy the regularization terms back into the Hessian
      for (size_t index: Range(number_regularization_nonzeros)) {
         hessian_values[number_hessian_nonzeros + index] = primal_regularization_values[index];
      }
   }

   void PrimalDualInertiaCorrection::regularize_hessian(Statistics& statistics, const Subproblem& subproblem,
         const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         View<double> primal_regularization_values) {
      this->primal_regularization = 0.;
      this->dual_regularization = 0.;
      size_t number_attempts = 0;

      while (true) {
         primal_regularization_values.fill(this->primal_regularization);
         DEBUG << "Testing factorization with primal regularization factor " << this->primal_regularization << '\n';
         linear_solver.do_numerical_factorization(subproblem.is_hessian_positive_definite());
         ++number_attempts;
         const Inertia estimated_inertia = linear_solver.get_inertia();
         DEBUG << "Expected inertia  " << expected_inertia << '\n';
         DEBUG << "Estimated inertia " << estimated_inertia << '\n';

         if (estimated_inertia == expected_inertia) {
            DEBUG << "The inertia is correct\n";
            if (0. < this->primal_regularization) {
               this->previous_primal_regularization = this->primal_regularization;
            }
            statistics.set("Prim reg", this->primal_regularization);
            statistics.set("Dual reg", this->dual_regularization);
            return;
         }
         this->increase_primal_regularization(number_attempts);
      }
   }

   void PrimalDualInertiaCorrection::regularize_augmented_matrix(Statistics& statistics, const Subproblem& /*subproblem*/,
         double dual_regularization_parameter, const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         BlockAugmentedMatrix block_augmented_matrix, Range<> primal_regularization_variables) {
      this->primal_regularization = 0.;
      this->dual_regularization = 0.;
      this->factorize_with_correct_inertia(statistics, dual_regularization_parameter, expected_inertia, linear_solver,
         block_augmented_matrix, primal_regularization_variables, /* assume_singular = */ false);
   }

   // the matrix was factorized with the current regularization factors, but is deemed singular by the caller (IPOPT's
   // "pretend singular" heuristic): resume the correction from the current factors as if the factorization were singular
   bool PrimalDualInertiaCorrection::regularize_singular_augmented_matrix(Statistics& statistics, const Subproblem& /*subproblem*/,
         double dual_regularization_parameter, const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         BlockAugmentedMatrix block_augmented_matrix, Range<> primal_regularization_variables) {
      DEBUG << "The augmented matrix is assumed singular\n";
      this->factorize_with_correct_inertia(statistics, dual_regularization_parameter, expected_inertia, linear_solver,
         block_augmented_matrix, primal_regularization_variables, /* assume_singular = */ true);
      return true;
   }

   // refactorize the matrix with the current regularization factors (e.g. after the pivoting tolerance of the linear solver
   // was increased): the inertia of the new factorization is checked and the regularization increased if needed
   void PrimalDualInertiaCorrection::refactorize_augmented_matrix(Statistics& statistics, const Subproblem& /*subproblem*/,
         double dual_regularization_parameter, const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         BlockAugmentedMatrix block_augmented_matrix, Range<> primal_regularization_variables) {
      this->factorize_with_correct_inertia(statistics, dual_regularization_parameter, expected_inertia, linear_solver,
         block_augmented_matrix, primal_regularization_variables, /* assume_singular = */ false);
   }

   void PrimalDualInertiaCorrection::factorize_with_correct_inertia(Statistics& statistics, double dual_regularization_parameter,
         const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         BlockAugmentedMatrix& block_augmented_matrix, Range<> primal_regularization_variables, bool assume_singular) {
      const double dual_regularization_value = this->dual_regularization_fraction * dual_regularization_parameter;
      size_t number_attempts = 0;

      while (true) {
         bool singular, missing_negative, missing_positive;
         if (assume_singular) {
            // the matrix was already factorized with the current factors: skip the factorization and treat it as singular
            assume_singular = false;
            singular = true;
            missing_negative = missing_positive = false;
         }
         else {
            if (this->regularize_all_variables) {
               block_augmented_matrix.primal_inertia_correction.fill(this->primal_regularization);
            }
            else {
               block_augmented_matrix.primal_inertia_correction.fill(0.);
               for (size_t variable_index: primal_regularization_variables) {
                  block_augmented_matrix.primal_inertia_correction[variable_index] = this->primal_regularization;
               }
            }
            block_augmented_matrix.dual_inertia_correction.fill(-this->dual_regularization);
            DEBUG << "Testing factorization with regularization factors (" << this->primal_regularization << ", " << this->dual_regularization << ")\n";
            linear_solver.do_numerical_factorization(false);
            ++number_attempts;
            const Inertia estimated_inertia = linear_solver.get_inertia();
            DEBUG << "Expected inertia  " << expected_inertia << '\n';
            DEBUG << "Estimated inertia " << estimated_inertia << '\n';

            if (estimated_inertia == expected_inertia) {
               DEBUG << "The inertia is correct\n";
               if (0. < this->primal_regularization) {
                  this->previous_primal_regularization = this->primal_regularization;
               }
               statistics.set("Prim reg", this->primal_regularization);
               statistics.set("Dual reg", this->dual_regularization);
               return;
            }
            singular = linear_solver.matrix_is_singular() || 0 < estimated_inertia.zero;
            missing_negative = estimated_inertia.negative < expected_inertia.negative;
            missing_positive = estimated_inertia.positive < expected_inertia.positive;
         }
         DEBUG << "Number of attempts: " << number_attempts << '\n';

         // missing negative eigenvalues (incl. zero ones or a surplus of positive ones): rank-deficient Jacobian -> dual regularization.
         // Primal regularization only adds positive eigenvalues and cannot fix this
         bool dual_regularization_added = false;
         if ((singular || missing_negative) && this->dual_regularization == 0. && 0. < dual_regularization_value) {
            DEBUG << "Adding dual regularization\n";
            this->dual_regularization = dual_regularization_value;
            dual_regularization_added = true;
         }

         // missing positive eigenvalues, or dual regularization is unavailable/exhausted -> increase the primal regularization
         if (missing_positive || !dual_regularization_added) {
            this->increase_primal_regularization(number_attempts);
         }
      }
   }

   // IPOPT-like update of the primal regularization after an unsuccessful factorization
   void PrimalDualInertiaCorrection::increase_primal_regularization(size_t number_attempts) {
      if (this->primal_regularization == 0.) {
         this->primal_regularization = (this->previous_primal_regularization == 0.) ? this->primal_regularization_initial_factor :
            std::max(this->primal_regularization_lb, this->previous_primal_regularization / this->primal_regularization_decrease_factor);
      }
      else if (this->previous_primal_regularization == 0. || this->threshold_unsuccessful_attempts < number_attempts) {
         this->primal_regularization *= this->primal_regularization_fast_increase_factor;
      }
      else {
         this->primal_regularization *= this->primal_regularization_slow_increase_factor;
      }
      if (this->regularization_failure_threshold < this->primal_regularization) {
         DEBUG << "The inertia correction failed\n";
         throw UnstableInertiaCorrection();
      }
   }

   bool PrimalDualInertiaCorrection::performs_primal_regularization() const {
      return true;
   }

   bool PrimalDualInertiaCorrection::performs_dual_regularization() const {
      return true;
   }

   [[nodiscard]] double PrimalDualInertiaCorrection::get_primal_regularization_factor() const {
      return this->primal_regularization;
   }

   std::string PrimalDualInertiaCorrection::get_name() const {
      return "primal-dual";
   }
} // namespace