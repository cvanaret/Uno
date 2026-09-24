// Copyright (c) 2025 Charlie Vanaret
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
         threshold_unsuccessful_attempts(options.get_unsigned_int("threshold_unsuccessful_attempts")) {
   }

   void PrimalDualInertiaCorrection::initialize_statistics(Statistics& statistics) {
      statistics.add_column("Prim reg", Statistics::double_width + 1, 2, /* is_extended = */ true);
      statistics.add_column("Dual reg", Statistics::double_width + 1, 2, /* is_extended = */ true);
   }

   void PrimalDualInertiaCorrection::regularize_hessian(Statistics& statistics, const Subproblem& subproblem,
         const Inertia& expected_inertia, View<double> hessian_values) {
      // pick the member linear solver
      if (this->optional_linear_solver == nullptr) {
         this->optional_linear_solver = SymmetricIndefiniteLinearSolverFactory::create(this->options);
         this->optional_linear_solver->get_linear_system().initialize_augmented_system(subproblem);
         this->optional_linear_solver->initialize_memory();
         this->optional_linear_solver->do_symbolic_analysis();
      }
      this->regularize_hessian(statistics, subproblem, expected_inertia, *this->optional_linear_solver, hessian_values);
   }

   void PrimalDualInertiaCorrection::regularize_hessian(Statistics& /*statistics*/, const Subproblem& /*subproblem*/,
         const Inertia& /*expected_inertia*/, DirectSymmetricIndefiniteLinearSolver<double>& /*linear_solver*/,
         View<double> /*hessian_values*/) {
      // to regularize the Hessian only, call the function for the augmented matrix with no dual part
      // TODO fix
      throw std::runtime_error("PrimalDualInertiaCorrection::regularize_hessian not implemented yet");
   }

   // the augmented matrix has been factorized prior to calling this function
   void PrimalDualInertiaCorrection::regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         double dual_regularization_parameter, const Inertia& expected_inertia, View<double> primal_inertia_correction_block,
         View<double> dual_inertia_correction_block) {
      if (this->optional_linear_solver == nullptr) {
         this->optional_linear_solver = SymmetricIndefiniteLinearSolverFactory::create(this->options);
         this->optional_linear_solver->get_linear_system().initialize_augmented_system(subproblem);
         this->optional_linear_solver->initialize_memory();
         this->optional_linear_solver->do_symbolic_analysis();
      }
      this->regularize_augmented_matrix(statistics, subproblem, dual_regularization_parameter, expected_inertia,
         *this->optional_linear_solver, primal_inertia_correction_block, dual_inertia_correction_block);
   }

   void PrimalDualInertiaCorrection::regularize_augmented_matrix(Statistics& statistics, const Subproblem& /*subproblem*/,
         double dual_regularization_parameter, const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         View<double> primal_inertia_correction_block, View<double> dual_inertia_correction_block) {
      const double dual_regularization_value = this->dual_regularization_fraction * dual_regularization_parameter;
      this->primal_regularization = 0.;
      this->dual_regularization = 0.;
      size_t number_attempts = 0;

      while (true) {
         primal_inertia_correction_block.fill(this->primal_regularization);
         dual_inertia_correction_block.fill(-this->dual_regularization);
         DEBUG << "Testing factorization with regularization factors (" << this->primal_regularization << ", " << this->dual_regularization << ")\n";
         linear_solver.do_numerical_factorization(false);
         ++number_attempts;
         const Inertia estimated_inertia = linear_solver.get_inertia();
         DEBUG << "Number of attempts: " << number_attempts << '\n';
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

         // missing negative eigenvalues (incl. zero ones or a surplus of positive ones): rank-deficient Jacobian -> dual regularization.
         // Primal regularization only adds positive eigenvalues and cannot fix this
         const bool singular = linear_solver.matrix_is_singular() || 0 < estimated_inertia.zero;
         const bool missing_negative = estimated_inertia.negative < expected_inertia.negative;
         bool dual_regularization_added = false;
         if ((singular || missing_negative) && this->dual_regularization == 0. && 0. < dual_regularization_value) {
            DEBUG << "Adding dual regularization\n";
            this->dual_regularization = dual_regularization_value;
            dual_regularization_added = true;
         }

         // missing positive eigenvalues, or dual regularization is unavailable/exhausted -> increase the primal regularization
         const bool missing_positive = estimated_inertia.positive < expected_inertia.positive;
         if (missing_positive || !dual_regularization_added) {
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