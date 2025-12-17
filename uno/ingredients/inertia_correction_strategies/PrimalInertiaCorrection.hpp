// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRIMALINERTIACORRECTION_H
#define UNO_PRIMALINERTIACORRECTION_H

#include <cassert>
#include <memory>
#include <string>
#include "InertiaCorrectionStrategy.hpp"
#include "UnstableRegularization.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   template <typename ElementType>
   class PrimalInertiaCorrection: public InertiaCorrectionStrategy<ElementType> {
   public:
      explicit PrimalInertiaCorrection(const Options& options);

      void initialize_statistics(Statistics& statistics) override;

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
         double* primal_regularization_values,
         double* dual_regularization_values) override;

      [[nodiscard]] bool performs_primal_regularization() const override;
      [[nodiscard]] bool performs_dual_regularization() const override;
      [[nodiscard]] double get_primal_regularization_factor() const override;
      [[nodiscard]] std::string get_name() const override;

   protected:
      const std::string& optional_linear_solver_name;
      std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<double>> optional_linear_solver{};
      double regularization_factor{0.};
      const double regularization_initial_value{};
      const double regularization_increase_factor{};
      const double regularization_failure_threshold{};
   };

   template <typename ElementType>
   PrimalInertiaCorrection<ElementType>::PrimalInertiaCorrection(const Options& options):
         InertiaCorrectionStrategy<ElementType>(),
         optional_linear_solver_name(options.get_string("linear_solver")),
         regularization_initial_value(options.get_double("regularization_initial_value")),
         regularization_increase_factor(options.get_double("regularization_increase_factor")),
         regularization_failure_threshold(options.get_double("regularization_failure_threshold")) {
   }

   template <typename ElementType>
   void PrimalInertiaCorrection<ElementType>::initialize_statistics(Statistics& statistics) {
      statistics.add_column("Regulariz", Statistics::double_width, 2, Statistics::column_order.at("Regulariz"));
   }

   // Nocedal and Wright, p51
   template <typename ElementType>
   void PrimalInertiaCorrection<ElementType>::regularize_hessian(Statistics& statistics, const Subproblem& subproblem,
         const double* hessian_values, const Inertia& expected_inertia, double* primal_regularization_values) {
      // pick the member linear solver
      if (this->optional_linear_solver == nullptr) {
         this->optional_linear_solver = SymmetricIndefiniteLinearSolverFactory::create(this->optional_linear_solver_name);
         this->optional_linear_solver->initialize_hessian(subproblem);
         this->optional_linear_solver->do_symbolic_analysis();
      }
      this->regularize_hessian(statistics, subproblem, hessian_values, expected_inertia, *this->optional_linear_solver,
         primal_regularization_values);
   }

   template <typename ElementType>
   void PrimalInertiaCorrection<ElementType>::regularize_hessian(Statistics& statistics, const Subproblem& subproblem,
         const double* hessian_values, const Inertia& expected_inertia,
         DirectSymmetricIndefiniteLinearSolver<double>& linear_solver, double* primal_regularization_values) {
      assert(hessian_values != nullptr);
      const size_t number_hessian_nonzeros = subproblem.number_regularized_hessian_nonzeros();

      this->regularization_factor = 0.;
      bool good_inertia = false;
      while (!good_inertia) {
         DEBUG << "Testing factorization with regularization factor " << this->regularization_factor << '\n';
         for (size_t index: Range(subproblem.get_primal_regularization_variables().size())) {
            primal_regularization_values[index] = this->regularization_factor;
         }
         DEBUG << "Current Hessian:";
         for (size_t nonzero_index: Range(number_hessian_nonzeros)) {
            DEBUG << ' ' << hessian_values[nonzero_index];
         }
         DEBUG << '\n';

         // perform factorization to get an estimate of the inertia
         linear_solver.do_numerical_factorization(hessian_values);

         // check inertia
         const Inertia estimated_inertia = linear_solver.get_inertia();
         DEBUG << "Expected inertia: " << expected_inertia << '\n';
         DEBUG << "Estimated inertia: " << estimated_inertia << '\n';
         if (estimated_inertia == expected_inertia) {
            good_inertia = true;
            DEBUG << "Factorization was a success\n";
         }
         else {
            this->regularization_factor = (this->regularization_factor == 0.) ? this->regularization_initial_value :
               this->regularization_increase_factor * this->regularization_factor;
            if (this->regularization_factor > this->regularization_failure_threshold) {
               throw UnstableRegularization();
            }
         }
      }
      statistics.set("Regulariz", this->regularization_factor);
   }

   template <typename ElementType>
   void PrimalInertiaCorrection<ElementType>::regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         const double* augmented_matrix_values, ElementType dual_regularization_parameter,
         const Inertia& expected_inertia, double* primal_regularization_values,
         double* dual_regularization_values) {
      // pick the member linear solver
      if (this->optional_linear_solver == nullptr) {
         this->optional_linear_solver = SymmetricIndefiniteLinearSolverFactory::create(this->optional_linear_solver_name);
         this->optional_linear_solver->initialize_hessian(subproblem);
         this->optional_linear_solver->do_symbolic_analysis();
      }
      this->regularize_augmented_matrix(statistics, subproblem, augmented_matrix_values, dual_regularization_parameter,
         expected_inertia, *this->optional_linear_solver, primal_regularization_values, dual_regularization_values);
   }

   template <typename ElementType>
   void PrimalInertiaCorrection<ElementType>::regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         const double* augmented_matrix_values, ElementType /*dual_regularization_parameter*/,
         const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         double* primal_regularization_values, double* /*dual_regularization_values*/) {
      this->regularize_hessian(statistics, subproblem, augmented_matrix_values, expected_inertia, linear_solver,
         primal_regularization_values);
   }

   template <typename ElementType>
   bool PrimalInertiaCorrection<ElementType>::performs_primal_regularization() const {
      return true;
   }

   template <typename ElementType>
   bool PrimalInertiaCorrection<ElementType>::performs_dual_regularization() const {
      return false;
   }

   template <typename ElementType>
   double PrimalInertiaCorrection<ElementType>::get_primal_regularization_factor() const {
      return this->regularization_factor;
   }

   template <typename ElementType>
   std::string PrimalInertiaCorrection<ElementType>::get_name() const {
      return "primal";
   }
} // namespace

#endif // UNO_PRIMALINERTIACORRECTION_H