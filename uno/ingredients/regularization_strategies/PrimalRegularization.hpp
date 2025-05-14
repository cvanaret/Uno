// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRIMALREGULARIZATION_H
#define UNO_PRIMALREGULARIZATION_H

#include <cassert>
#include <memory>
#include <string>
#include "RegularizationStrategy.hpp"
#include "UnstableRegularization.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   template <typename ElementType>
   class PrimalRegularization: public RegularizationStrategy<ElementType> {
   public:
      explicit PrimalRegularization(const Options& options);

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
   PrimalRegularization<ElementType>::PrimalRegularization(const Options& options):
         RegularizationStrategy<ElementType>(),
         optional_linear_solver_name(options.get_string("linear_solver")),
         regularization_initial_value(options.get_double("regularization_initial_value")),
         regularization_increase_factor(options.get_double("regularization_increase_factor")),
         regularization_failure_threshold(options.get_double("regularization_failure_threshold")) {
   }

   template <typename ElementType>
   void PrimalRegularization<ElementType>::initialize_statistics(Statistics& statistics, const Options& options) {
      statistics.add_column("regulariz", Statistics::double_width - 4, options.get_int("statistics_regularization_column_order"));
   }

   // Nocedal and Wright, p51
   template <typename ElementType>
   void PrimalRegularization<ElementType>::regularize_hessian(Statistics& statistics, const Subproblem& subproblem,
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
   void PrimalRegularization<ElementType>::regularize_hessian(Statistics& statistics, const Subproblem& subproblem,
         const double* hessian_values, const Inertia& expected_inertia,
         DirectSymmetricIndefiniteLinearSolver<double>& linear_solver, double* primal_regularization_values) {
      assert(hessian_values != nullptr);

      DEBUG << "Current Hessian:\n" << hessian_values << '\n';
      const double smallest_diagonal_entry = 0.; // TODO hessian_values.smallest_diagonal_entry(expected_inertia.positive);
      DEBUG << "The minimal diagonal entry of the matrix is " << smallest_diagonal_entry << '\n';

      this->regularization_factor = (smallest_diagonal_entry > 0.) ? 0. : this->regularization_initial_value - smallest_diagonal_entry;
      bool good_inertia = false;
      while (!good_inertia) {
         DEBUG << "Testing factorization with regularization factor " << this->regularization_factor << '\n';
         for (size_t index: Range(subproblem.get_primal_regularization_variables().size())) {
            primal_regularization_values[index] = this->regularization_factor;
         }
         DEBUG << "Current Hessian:\n" << hessian_values;

         linear_solver.do_numerical_factorization(hessian_values);
         const Inertia estimated_inertia = linear_solver.get_inertia();
         DEBUG << "Expected inertia: " << expected_inertia << '\n';
         DEBUG << "Estimated inertia: " << estimated_inertia << '\n';
         if (estimated_inertia == expected_inertia) {
            good_inertia = true;
            DEBUG << "Factorization was a success";
         }
         else {
            this->regularization_factor = (this->regularization_factor == 0.) ? this->regularization_initial_value :
               this->regularization_increase_factor * this->regularization_factor;
            if (this->regularization_factor > this->regularization_failure_threshold) {
               throw UnstableRegularization();
            }
         }
         DEBUG << '\n';
      }
      statistics.set("regulariz", this->regularization_factor);
   }

   template <typename ElementType>
   void PrimalRegularization<ElementType>::regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
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
   void PrimalRegularization<ElementType>::regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         const double* augmented_matrix_values, ElementType /*dual_regularization_parameter*/,
         const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         double* primal_regularization_values, double* /*dual_regularization_values*/) {
      this->regularize_hessian(statistics, subproblem, augmented_matrix_values, expected_inertia, linear_solver,
         primal_regularization_values);
   }

   template <typename ElementType>
   bool PrimalRegularization<ElementType>::performs_primal_regularization() const {
      return true;
   }

   template <typename ElementType>
   bool PrimalRegularization<ElementType>::performs_dual_regularization() const {
      return false;
   }

   template <typename ElementType>
   double PrimalRegularization<ElementType>::get_primal_regularization_factor() const {
      return this->regularization_factor;
   }

   template <typename ElementType>
   std::string PrimalRegularization<ElementType>::get_name() const {
      return "primal";
   }
} // namespace

#endif // UNO_PRIMALREGULARIZATION_H