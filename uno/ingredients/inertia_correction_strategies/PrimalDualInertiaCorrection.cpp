// Copyright (c) 2025-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
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

      while (true) {
         primal_regularization_values.fill(this->primal_regularization);
         DEBUG << "Testing factorization with primal regularization factor " << this->primal_regularization << '\n';
         linear_solver.do_numerical_factorization(subproblem.is_hessian_positive_definite());
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
         if (!this->increase_primal_regularization()) {
            DEBUG << "The inertia correction failed\n";
            throw UnstableInertiaCorrection();
         }
      }
   }

   // new matrix (IPOPT's ConsiderNewSystem): the regularization starts from 0, except for the blocks found to be
   // structurally degenerate on the previous matrices
   void PrimalDualInertiaCorrection::regularize_augmented_matrix(Statistics& statistics, const Subproblem& /*subproblem*/,
         double dual_regularization_parameter, const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         BlockAugmentedMatrix block_augmented_matrix, Range<> primal_regularization_variables) {
      const double dual_regularization_value = this->dual_regularization_fraction * dual_regularization_parameter;
      this->finalize_degeneracy_test();
      // remember the last nonzero primal regularization
      if (0. < this->primal_regularization) {
         this->previous_primal_regularization = this->primal_regularization;
      }
      const bool degeneracy_undetermined = (this->hessian_degeneracy == Degeneracy::NOT_YET_DETERMINED ||
         this->jacobian_degeneracy == Degeneracy::NOT_YET_DETERMINED);
      this->degeneracy_test = degeneracy_undetermined ? DegeneracyTest::NO_REGULARIZATION : DegeneracyTest::NONE;
      this->dual_regularization = (this->jacobian_degeneracy == Degeneracy::DEGENERATE) ? dual_regularization_value : 0.;
      this->primal_regularization = 0.;
      if (this->hessian_degeneracy == Degeneracy::DEGENERATE && !this->increase_primal_regularization()) {
         DEBUG << "The inertia correction failed\n";
         throw UnstableInertiaCorrection();
      }
      this->quality_improved_for_current_matrix = false;
      this->regularize_until_correct_inertia(statistics, dual_regularization_value, expected_inertia, linear_solver,
         block_augmented_matrix, primal_regularization_variables, /* assume_singular = */ false);
   }

   // the matrix was factorized with the current regularization factors, but is deemed singular by the caller (IPOPT's
   // "pretend singular" heuristic): resume the correction from the current factors as if the factorization were singular
   bool PrimalDualInertiaCorrection::regularize_singular_augmented_matrix(Statistics& statistics, const Subproblem& /*subproblem*/,
         double dual_regularization_parameter, const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         BlockAugmentedMatrix block_augmented_matrix, Range<> primal_regularization_variables) {
      DEBUG << "The augmented matrix is assumed singular\n";
      this->regularize_until_correct_inertia(statistics, this->dual_regularization_fraction * dual_regularization_parameter,
         expected_inertia, linear_solver, block_augmented_matrix, primal_regularization_variables, /* assume_singular = */ true);
      return true;
   }

   // refactorize the matrix with the current regularization factors (e.g. after the pivoting tolerance of the linear solver
   // was increased): the inertia of the new factorization is checked and the regularization increased if needed
   void PrimalDualInertiaCorrection::refactorize_augmented_matrix(Statistics& statistics, const Subproblem& /*subproblem*/,
         double dual_regularization_parameter, const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         BlockAugmentedMatrix block_augmented_matrix, Range<> primal_regularization_variables) {
      this->regularize_until_correct_inertia(statistics, this->dual_regularization_fraction * dual_regularization_parameter,
         expected_inertia, linear_solver, block_augmented_matrix, primal_regularization_variables, /* assume_singular = */ false);
   }

   // IPOPT's PDFullSpaceSolver::SolveOnce
   void PrimalDualInertiaCorrection::regularize_until_correct_inertia(Statistics& statistics, double dual_regularization_value,
         const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         BlockAugmentedMatrix& block_augmented_matrix, Range<> primal_regularization_variables, bool assume_singular) {
      const bool has_constraints = (0 < expected_inertia.negative);
      while (true) {
         if (assume_singular) {
            // the matrix was already factorized with the current factors: skip the factorization and treat it as singular
            assume_singular = false;
            this->regularize_for_singularity(dual_regularization_value);
            continue;
         }
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
         const Inertia estimated_inertia = linear_solver.get_inertia();
         DEBUG << "Expected inertia  " << expected_inertia << '\n';
         DEBUG << "Estimated inertia " << estimated_inertia << '\n';

         if (estimated_inertia == expected_inertia) {
            DEBUG << "The inertia is correct\n";
            statistics.set("Prim reg", this->primal_regularization);
            statistics.set("Dual reg", this->dual_regularization);
            return;
         }
         const bool singular = linear_solver.matrix_is_singular() || 0 < estimated_inertia.zero;
         if (singular && has_constraints) {
            this->regularize_for_singularity(dual_regularization_value);
         }
         else if (!singular && estimated_inertia.negative < expected_inertia.negative) {
            // too few negative eigenvalues: improve the quality of the factorization once, then assume the matrix is singular
            if (!this->quality_improved_for_current_matrix && linear_solver.increase_pivoting_tolerance()) {
               DEBUG << "Too few negative eigenvalues: refactorizing with an increased pivoting tolerance\n";
               this->quality_improved_for_current_matrix = true;
            }
            else {
               DEBUG << "Too few negative eigenvalues: the matrix is assumed singular\n";
               this->regularize_for_singularity(dual_regularization_value);
            }
         }
         else {
            this->regularize_for_wrong_inertia(dual_regularization_value);
         }
      }
   }

   // IPOPT's PerturbForSingularity: while the degeneracy of the Hessian and the Jacobian is not determined, the
   // dual regularization alone, then the primal regularization alone, then both are tested
   void PrimalDualInertiaCorrection::regularize_for_singularity(double dual_regularization_value) {
      DEBUG << "Regularizing a singular matrix\n";
      bool success = true;
      if (this->hessian_degeneracy == Degeneracy::NOT_YET_DETERMINED || this->jacobian_degeneracy == Degeneracy::NOT_YET_DETERMINED) {
         switch (this->degeneracy_test) {
            case DegeneracyTest::NO_REGULARIZATION:
               if (this->jacobian_degeneracy == Degeneracy::NOT_YET_DETERMINED) {
                  this->dual_regularization = dual_regularization_value;
                  this->degeneracy_test = DegeneracyTest::DUAL_ONLY;
               }
               else {
                  success = this->increase_primal_regularization();
                  this->degeneracy_test = DegeneracyTest::PRIMAL_ONLY;
               }
               break;
            case DegeneracyTest::DUAL_ONLY:
               this->dual_regularization = 0.;
               success = this->increase_primal_regularization();
               this->degeneracy_test = DegeneracyTest::PRIMAL_ONLY;
               break;
            case DegeneracyTest::PRIMAL_ONLY:
               this->dual_regularization = dual_regularization_value;
               success = this->increase_primal_regularization();
               this->degeneracy_test = DegeneracyTest::PRIMAL_AND_DUAL;
               break;
            case DegeneracyTest::PRIMAL_AND_DUAL:
            case DegeneracyTest::NONE:
               success = this->increase_primal_regularization();
               break;
         }
      }
      // the constraints were already regularized: treat the matrix as having the wrong inertia
      else if (0. < this->dual_regularization) {
         success = this->increase_primal_regularization();
      }
      else {
         this->dual_regularization = dual_regularization_value;
      }
      if (!success) {
         DEBUG << "The inertia correction failed\n";
         throw UnstableInertiaCorrection();
      }
   }

   // IPOPT's PerturbForWrongInertia: increase the primal regularization. If it becomes too large without dual
   // regularization, the primal regularization starts over with dual regularization
   void PrimalDualInertiaCorrection::regularize_for_wrong_inertia(double dual_regularization_value) {
      this->finalize_degeneracy_test();
      bool success = this->increase_primal_regularization();
      if (!success && this->dual_regularization == 0.) {
         this->dual_regularization = dual_regularization_value;
         this->primal_regularization = 0.;
         this->degeneracy_test = DegeneracyTest::NONE;
         if (this->hessian_degeneracy == Degeneracy::DEGENERATE) {
            this->hessian_degeneracy = Degeneracy::NOT_YET_DETERMINED;
         }
         success = this->increase_primal_regularization();
      }
      if (!success) {
         DEBUG << "The inertia correction failed\n";
         throw UnstableInertiaCorrection();
      }
   }

   // IPOPT's get_deltas_for_wrong_inertia: the first nonzero primal regularization is a fraction of the last nonzero one,
   // the subsequent ones increase fast if there is no previous one or if the regularization is far above it, slowly
   // otherwise. Returns false if the regularization exceeds its maximum
   bool PrimalDualInertiaCorrection::increase_primal_regularization() {
      if (this->primal_regularization == 0.) {
         this->primal_regularization = (this->previous_primal_regularization == 0.) ? this->primal_regularization_initial_factor :
            std::max(this->primal_regularization_lb, this->previous_primal_regularization / this->primal_regularization_decrease_factor);
      }
      else if (this->previous_primal_regularization == 0. || 1e5 * this->previous_primal_regularization < this->primal_regularization) {
         this->primal_regularization *= this->primal_regularization_fast_increase_factor;
      }
      else {
         this->primal_regularization *= this->primal_regularization_slow_increase_factor;
      }
      if (this->regularization_failure_threshold < this->primal_regularization) {
         DEBUG << "The primal regularization is too large: " << this->primal_regularization << '\n';
         this->previous_primal_regularization = 0.;
         return false;
      }
      return true;
   }

   // IPOPT's finalize_test: conclude from the regularization that made the previous factorization nonsingular. A block
   // is deemed degenerate after it required regularization on degeneracy_threshold matrices
   void PrimalDualInertiaCorrection::finalize_degeneracy_test() {
      switch (this->degeneracy_test) {
         case DegeneracyTest::NONE:
            return;
         case DegeneracyTest::NO_REGULARIZATION:
            if (this->hessian_degeneracy == Degeneracy::NOT_YET_DETERMINED) {
               this->hessian_degeneracy = Degeneracy::NOT_DEGENERATE;
            }
            if (this->jacobian_degeneracy == Degeneracy::NOT_YET_DETERMINED) {
               this->jacobian_degeneracy = Degeneracy::NOT_DEGENERATE;
            }
            break;
         case DegeneracyTest::DUAL_ONLY:
            if (this->hessian_degeneracy == Degeneracy::NOT_YET_DETERMINED) {
               this->hessian_degeneracy = Degeneracy::NOT_DEGENERATE;
            }
            if (this->jacobian_degeneracy == Degeneracy::NOT_YET_DETERMINED) {
               ++this->number_degenerate_matrices;
               if (this->degeneracy_threshold <= this->number_degenerate_matrices) {
                  DEBUG << "The Jacobian is deemed degenerate\n";
                  this->jacobian_degeneracy = Degeneracy::DEGENERATE;
               }
            }
            break;
         case DegeneracyTest::PRIMAL_ONLY:
            if (this->jacobian_degeneracy == Degeneracy::NOT_YET_DETERMINED) {
               this->jacobian_degeneracy = Degeneracy::NOT_DEGENERATE;
            }
            if (this->hessian_degeneracy == Degeneracy::NOT_YET_DETERMINED) {
               ++this->number_degenerate_matrices;
               if (this->degeneracy_threshold <= this->number_degenerate_matrices) {
                  DEBUG << "The Hessian is deemed degenerate\n";
                  this->hessian_degeneracy = Degeneracy::DEGENERATE;
               }
            }
            break;
         case DegeneracyTest::PRIMAL_AND_DUAL:
            ++this->number_degenerate_matrices;
            if (this->degeneracy_threshold <= this->number_degenerate_matrices) {
               DEBUG << "The Hessian and the Jacobian are deemed degenerate\n";
               this->hessian_degeneracy = Degeneracy::DEGENERATE;
               this->jacobian_degeneracy = Degeneracy::DEGENERATE;
            }
            break;
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