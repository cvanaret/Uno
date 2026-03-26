// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "PrimalInertiaCorrection.hpp"
#include "UnstableInertiaCorrection.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem_solvers/LinearSystem.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   PrimalInertiaCorrection::PrimalInertiaCorrection(const Options& options):
         InertiaCorrectionStrategy(),
         optional_linear_solver_name(options.get_string("linear_solver")),
         regularization_initial_factor(options.get_double("primal_regularization_initial_factor")),
         regularization_increase_factor(options.get_double("regularization_increase_factor")),
         regularization_failure_threshold(options.get_double("regularization_failure_threshold")) {
   }

   void PrimalInertiaCorrection::initialize_statistics(Statistics& statistics) {
      statistics.add_column("Regulariz", Statistics::double_width, 2, Statistics::column_order.at("Regulariz"));
   }

   // Nocedal and Wright, p51
   void PrimalInertiaCorrection::regularize_hessian(Statistics& statistics, const Subproblem& subproblem,
         const Inertia& expected_inertia, double* hessian_values) {
      // pick the member linear solver
      if (this->optional_linear_solver == nullptr) {
         this->optional_linear_solver = SymmetricIndefiniteLinearSolverFactory::create(this->optional_linear_solver_name);
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
      // figure out where to set the regularization terms in the linear system
      double* primal_regularization_values = matrix + number_hessian_nonzeros;
      // regularize the Hessian
      this->regularize_hessian(statistics, subproblem, expected_inertia, *this->optional_linear_solver, primal_regularization_values);
      // copy the regularization terms back into the Hessian
      for (size_t index: Range(subproblem.get_primal_regularization_variables().size())) {
         hessian_values[number_hessian_nonzeros + index] = primal_regularization_values[index];
      }
   }

   void PrimalInertiaCorrection::regularize_hessian(Statistics& statistics, const Subproblem& subproblem,
         const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         double* primal_regularization_values) {
      this->regularization_factor = 0.;
      bool good_inertia = false;
      while (!good_inertia) {
         DEBUG << "Testing factorization with regularization factor " << this->regularization_factor << '\n';
         for (size_t index: Range(subproblem.get_primal_regularization_variables().size())) {
            primal_regularization_values[index] = this->regularization_factor;
         }
         DEBUG << '\n';

         // perform factorization to get an estimate of the inertia
         const bool is_matrix_positive_definite = subproblem.is_hessian_positive_definite();
         linear_solver.do_numerical_factorization(is_matrix_positive_definite);

         // check inertia
         const Inertia estimated_inertia = linear_solver.get_inertia();
         DEBUG << "Expected inertia: " << expected_inertia << '\n';
         DEBUG << "Estimated inertia: " << estimated_inertia << '\n';
         if (estimated_inertia == expected_inertia) {
            good_inertia = true;
            DEBUG << "Factorization was a success\n";
         }
         else {
            this->regularization_factor = (this->regularization_factor == 0.) ? this->regularization_initial_factor :
               this->regularization_increase_factor * this->regularization_factor;
            if (this->regularization_factor > this->regularization_failure_threshold) {
               throw UnstableInertiaCorrection();
            }
         }
      }
      statistics.set("Regulariz", this->regularization_factor);
   }

   void PrimalInertiaCorrection::regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         double dual_regularization_parameter, const Inertia& expected_inertia, double* primal_regularization_values,
         double* dual_regularization_values) {
      // pick the member linear solver
      if (this->optional_linear_solver == nullptr) {
         this->optional_linear_solver = SymmetricIndefiniteLinearSolverFactory::create(this->optional_linear_solver_name);
         this->optional_linear_solver->get_linear_system().initialize_hessian(subproblem);
         this->optional_linear_solver->initialize_memory();
         this->optional_linear_solver->do_symbolic_analysis();
      }
      this->regularize_augmented_matrix(statistics, subproblem, dual_regularization_parameter, expected_inertia,
         *this->optional_linear_solver, primal_regularization_values, dual_regularization_values);
   }

   void PrimalInertiaCorrection::regularize_augmented_matrix(Statistics& statistics, const Subproblem& subproblem,
         double /*dual_regularization_parameter*/, const Inertia& expected_inertia, DirectSymmetricIndefiniteLinearSolver<double>& linear_solver,
         double* primal_regularization_values, double* /*dual_regularization_values*/) {
      this->regularize_hessian(statistics, subproblem, expected_inertia, linear_solver, primal_regularization_values);
   }

   bool PrimalInertiaCorrection::performs_primal_regularization() const {
      return true;
   }

   bool PrimalInertiaCorrection::performs_dual_regularization() const {
      return false;
   }

   double PrimalInertiaCorrection::get_primal_regularization_factor() const {
      return this->regularization_factor;
   }

   std::string PrimalInertiaCorrection::get_name() const {
      return "primal";
   }
} // namespace