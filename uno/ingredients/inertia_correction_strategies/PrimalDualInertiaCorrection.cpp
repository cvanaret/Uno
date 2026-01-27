// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "PrimalDualInertiaCorrection.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "options/Options.hpp"
#include "symbolic/Collection.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   PrimalDualInertiaCorrection::PrimalDualInertiaCorrection(const Options& options):
         InertiaCorrectionStrategy(),
         optional_linear_solver_name(options.get_string("linear_solver")),
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
      statistics.add_column("Regulariz", Statistics::double_width, 2, Statistics::column_order.at("Regulariz"));
   }

   void PrimalDualInertiaCorrection::regularize_hessian(Statistics& statistics, const Subproblem& subproblem,
         const double* hessian_values, const Inertia& expected_inertia, double* primal_regularization_values) {
      // pick the member linear solver
      if (this->optional_linear_solver == nullptr) {
         this->optional_linear_solver = SymmetricIndefiniteLinearSolverFactory::create(this->optional_linear_solver_name);
         this->optional_linear_solver->initialize_augmented_system(subproblem);
         this->optional_linear_solver->do_symbolic_analysis();
      }
      this->regularize_hessian(statistics, subproblem, hessian_values, expected_inertia, *this->optional_linear_solver,
         primal_regularization_values);
   }

   void PrimalDualInertiaCorrection::regularize_hessian(Statistics& /*statistics*/, const Subproblem& /*subproblem*/,
         const double* /*hessian_values*/, const Inertia& /*expected_inertia*/,
         DirectSymmetricIndefiniteLinearSolver<double>& /*linear_solver*/,
         double* /*primal_regularization_values*/) {
      // to regularize the Hessian only, call the function for the augmented matrix with no dual part
      // TODO fix
      throw std::runtime_error("PrimalDualInertiaCorrection::regularize_hessian not implemented yet");
   }

   // the augmented matrix has been factorized prior to calling this function
   void PrimalDualInertiaCorrection::regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         const double* augmented_matrix_values, double dual_regularization_parameter,
         const Inertia& expected_inertia, double* primal_regularization_values, double* dual_regularization_values) {
      if (this->optional_linear_solver == nullptr) {
         this->optional_linear_solver = SymmetricIndefiniteLinearSolverFactory::create(this->optional_linear_solver_name);
         this->optional_linear_solver->initialize_augmented_system(subproblem);
         this->optional_linear_solver->do_symbolic_analysis();
      }
      this->regularize_augmented_matrix(statistics, subproblem, augmented_matrix_values, dual_regularization_parameter,
         expected_inertia, *this->optional_linear_solver, primal_regularization_values, dual_regularization_values);
   }

   void PrimalDualInertiaCorrection::regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         const double* augmented_matrix_values, double dual_regularization_parameter,
         const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         double* primal_regularization_values, double* dual_regularization_values) {
      assert(augmented_matrix_values != nullptr);

      this->primal_regularization = 0.;
      this->dual_regularization = 0.;
      for (size_t index: Range(subproblem.get_primal_regularization_variables().size())) {
         primal_regularization_values[index] = this->primal_regularization;
      }
      for (size_t index: Range(subproblem.get_dual_regularization_constraints().size())) {
         dual_regularization_values[index] = -this->dual_regularization;
      }
      DEBUG2 << "Original matrix values\n" << augmented_matrix_values << '\n';
      DEBUG << "Testing factorization with regularization factors (0, 0)\n";
      size_t number_attempts = 1;
      DEBUG << "Number of attempts: " << number_attempts << "\n\n";

      DEBUG << "Performing numerical factorization of the indefinite system\n";
      linear_solver.do_numerical_factorization(augmented_matrix_values);
      const Inertia estimated_inertia = linear_solver.get_inertia();
      DEBUG << "Expected inertia  " << expected_inertia << '\n';
      DEBUG << "Estimated inertia " << estimated_inertia << '\n';

      if (estimated_inertia == expected_inertia) {
         DEBUG << "The inertia is correct\n";
         statistics.set("Regulariz", this->primal_regularization);
         return;
      }

      // set the constraint regularization coefficient
      if (linear_solver.matrix_is_singular()) {
         DEBUG << "Matrix is singular\n";
         this->dual_regularization = this->dual_regularization_fraction * dual_regularization_parameter;
      }
      // set the Hessian regularization coefficient
      if (this->previous_primal_regularization == 0.) {
         this->primal_regularization = this->primal_regularization_initial_factor;
      }
      else {
         this->primal_regularization = std::max(this->primal_regularization_lb,
            this->previous_primal_regularization / this->primal_regularization_decrease_factor);
      }

      // regularize the augmented matrix
      for (size_t index: Range(subproblem.get_primal_regularization_variables().size())) {
         primal_regularization_values[index] = this->primal_regularization;
      }
      for (size_t index: Range(subproblem.get_dual_regularization_constraints().size())) {
         dual_regularization_values[index] = -this->dual_regularization;
      }

      bool good_inertia = false;
      while (!good_inertia) {
         DEBUG << "Testing factorization with regularization factors (" << this->primal_regularization << ", " << this->dual_regularization << ")\n";
         DEBUG2 << augmented_matrix_values << '\n';
         DEBUG << "Performing numerical factorization of the indefinite system\n";
         linear_solver.do_numerical_factorization(augmented_matrix_values);
         ++number_attempts;
         DEBUG << "Number of attempts: " << number_attempts << "\n";

         const Inertia new_estimated_inertia = linear_solver.get_inertia();
         DEBUG << "Expected inertia  " << expected_inertia << '\n';
         DEBUG << "Estimated inertia " << new_estimated_inertia << '\n';

         if (new_estimated_inertia.positive == expected_inertia.positive && new_estimated_inertia.negative == expected_inertia.negative &&
               new_estimated_inertia.zero == expected_inertia.zero) {
            good_inertia = true;
            DEBUG << "The inertia is correct\n";
            this->previous_primal_regularization = this->primal_regularization;
         }
         else {
            if (this->previous_primal_regularization == 0. || this->threshold_unsuccessful_attempts < number_attempts) {
               this->primal_regularization *= this->primal_regularization_fast_increase_factor;
            }
            else {
               this->primal_regularization *= this->primal_regularization_slow_increase_factor;
            }

            if (this->primal_regularization <= this->regularization_failure_threshold) {
               // regularize the augmented matrix
               for (size_t index: Range(subproblem.get_primal_regularization_variables().size())) {
                  primal_regularization_values[index] = this->primal_regularization;
               }
               for (size_t index: Range(subproblem.get_dual_regularization_constraints().size())) {
                  dual_regularization_values[index] = -this->dual_regularization;
               }
            }
            else {
               throw UnstableRegularization();
            }
         }
      }
      statistics.set("Regulariz", this->primal_regularization);
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