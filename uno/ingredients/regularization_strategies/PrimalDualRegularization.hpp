// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRIMALDUALREGULARIZATION_H
#define UNO_PRIMALDUALREGULARIZATION_H

#include <cassert>
#include <string>
#include "RegularizationStrategy.hpp"
#include "UnstableRegularization.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "options/Options.hpp"
#include "symbolic/Collection.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   template <typename ElementType>
   class PrimalDualRegularization: public RegularizationStrategy<ElementType> {
   public:
      explicit PrimalDualRegularization(const Options& options);

      void initialize_statistics(Statistics& statistics, const Options& options) override;

      void regularize_hessian(Statistics& statistics, const Subproblem& subproblem, const double* hessian_values,
         const Inertia& expected_inertia, double* primal_regularization_values) override;
      void regularize_hessian(Statistics& statistics, const Subproblem& subproblem, const double* hessian_values,
         const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         double* primal_regularization_values) override;
      void regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         const double* augmented_matrix_values, ElementType dual_regularization_parameter,
         const Inertia& expected_inertia, double* primal_regularization_values,
         double* dual_regularization_values) override;
      void regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         const double* augmented_matrix_values, ElementType dual_regularization_parameter,
         const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         double* primal_regularization_values, double* dual_regularization_values) override;

      [[nodiscard]] bool performs_primal_regularization() const override;
      [[nodiscard]] bool performs_dual_regularization() const override;
      [[nodiscard]] double get_primal_regularization_factor() const override;

      [[nodiscard]] std::string get_name() const override;

   protected:
      const std::string& optional_linear_solver_name;
      std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<double>> optional_linear_solver{};
      ElementType primal_regularization{0.};
      ElementType dual_regularization{0.};
      ElementType previous_primal_regularization{0.};
      const ElementType regularization_failure_threshold;
      const ElementType primal_regularization_initial_factor;
      const ElementType dual_regularization_fraction;
      const ElementType primal_regularization_lb;
      const ElementType primal_regularization_decrease_factor;
      const ElementType primal_regularization_fast_increase_factor;
      const ElementType primal_regularization_slow_increase_factor;
      const size_t threshold_unsuccessful_attempts;
   };

   template <typename ElementType>
   PrimalDualRegularization<ElementType>::PrimalDualRegularization(const Options& options):
         RegularizationStrategy<ElementType>(),
         optional_linear_solver_name(options.get_string("linear_solver")),
         regularization_failure_threshold(ElementType(options.get_double("regularization_failure_threshold"))),
         primal_regularization_initial_factor(ElementType(options.get_double("primal_regularization_initial_factor"))),
         dual_regularization_fraction(ElementType(options.get_double("dual_regularization_fraction"))),
         primal_regularization_lb(ElementType(options.get_double("primal_regularization_lb"))),
         primal_regularization_decrease_factor(ElementType(options.get_double("primal_regularization_decrease_factor"))),
         primal_regularization_fast_increase_factor(ElementType(options.get_double("primal_regularization_fast_increase_factor"))),
         primal_regularization_slow_increase_factor(ElementType(options.get_double("primal_regularization_slow_increase_factor"))),
         threshold_unsuccessful_attempts(options.get_unsigned_int("threshold_unsuccessful_attempts")) {
   }

   template <typename ElementType>
   void PrimalDualRegularization<ElementType>::initialize_statistics(Statistics& statistics, const Options& options) {
      statistics.add_column("regulariz", Statistics::double_width - 4, options.get_int("statistics_regularization_column_order"));
   }

   template <typename ElementType>
   void PrimalDualRegularization<ElementType>::regularize_hessian(Statistics& statistics, const Subproblem& subproblem,
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

   template <typename ElementType>
   void PrimalDualRegularization<ElementType>::regularize_hessian(Statistics& /*statistics*/, const Subproblem& /*subproblem*/,
         const double* /*hessian_values*/, const Inertia& /*expected_inertia*/,
         DirectSymmetricIndefiniteLinearSolver<double>& /*linear_solver*/,
         double* /*primal_regularization_values*/) {
      // to regularize the Hessian only, call the function for the augmented matrix with no dual part
      // TODO fix
      throw std::runtime_error("PrimalDualRegularization::regularize_hessian not implemented yet");
   }

   // the augmented matrix has been factorized prior to calling this function
   template <typename ElementType>
   void PrimalDualRegularization<ElementType>::regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         const double* augmented_matrix_values, ElementType dual_regularization_parameter,
         const Inertia& expected_inertia, double* primal_regularization_values, double* dual_regularization_values) {
      if (this->optional_linear_solver == nullptr) {
         this->optional_linear_solver = SymmetricIndefiniteLinearSolverFactory::create(this->optional_linear_solver_name);
         this->optional_linear_solver->initialize_augmented_system(subproblem);
         this->optional_linear_solver->do_symbolic_analysis();
      }
      this->regularize_augmented_matrix(statistics, subproblem, augmented_matrix_values, dual_regularization_parameter,
         expected_inertia, *this->optional_linear_solver, primal_regularization_values, dual_regularization_values);
   }

   template <typename ElementType>
   void PrimalDualRegularization<ElementType>::regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         const double* augmented_matrix_values, ElementType dual_regularization_parameter,
         const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         double* primal_regularization_values, double* dual_regularization_values) {
      assert(augmented_matrix_values != nullptr);

      this->primal_regularization = ElementType(0);
      this->dual_regularization = ElementType(0);
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
         statistics.set("regulariz", this->primal_regularization);
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
      statistics.set("regulariz", this->primal_regularization);
   }

   template <typename ElementType>
   bool PrimalDualRegularization<ElementType>::performs_primal_regularization() const {
      return true;
   }

   template <typename ElementType>
   bool PrimalDualRegularization<ElementType>::performs_dual_regularization() const {
      return true;
   }

   template <typename ElementType>
   [[nodiscard]] double PrimalDualRegularization<ElementType>::get_primal_regularization_factor() const {
      return this->primal_regularization;
   }

   template <typename ElementType>
   std::string PrimalDualRegularization<ElementType>::get_name() const {
      return "primal-dual";
   }
} // namespace

#endif // UNO_PRIMALDUALREGULARIZATION_H