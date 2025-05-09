// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRIMALREGULARIZATION_H
#define UNO_PRIMALREGULARIZATION_H

#include <memory>
#include "RegularizationStrategy.hpp"
#include "UnstableRegularization.hpp"
#include "ingredients/constraint_relaxation_strategies/OptimizationProblem.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   template <typename ElementType>
   class PrimalRegularization: public RegularizationStrategy<ElementType> {
   public:
      explicit PrimalRegularization(const Options& options);

      void initialize_memory(const OptimizationProblem& problem, const HessianModel& hessian_model) override;
      void initialize_statistics(Statistics& statistics, const Options& options) override;

      void regularize_hessian(Statistics& statistics, SymmetricMatrix<size_t, ElementType>& hessian, const Inertia& expected_inertia) override;
      void regularize_augmented_matrix(Statistics& statistics, SymmetricMatrix<size_t, ElementType>& augmented_matrix,
         const Collection<size_t>& primal_variables, const Collection<size_t>& dual_variables, ElementType dual_regularization_parameter,
         const Inertia& expected_inertia) override;

   protected:
      std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<size_t, double>> linear_solver;
      const double regularization_initial_value{};
      const double regularization_increase_factor{};
      const double regularization_failure_threshold{};
   };

   template <typename ElementType>
   PrimalRegularization<ElementType>::PrimalRegularization(const Options& options):
         RegularizationStrategy<ElementType>(),
         linear_solver(SymmetricIndefiniteLinearSolverFactory::create(options)),
         regularization_initial_value(options.get_double("regularization_initial_value")),
         regularization_increase_factor(options.get_double("regularization_increase_factor")),
         regularization_failure_threshold(options.get_double("regularization_failure_threshold")) {
   }

   template <typename ElementType>
   void PrimalRegularization<ElementType>::initialize_memory(const OptimizationProblem& problem, const HessianModel& hessian_model) {
      const size_t number_hessian_nonzeros = problem.number_hessian_nonzeros(hessian_model) +
         problem.number_variables; // diagonal regularization
      this->linear_solver->initialize_memory(problem.number_variables, number_hessian_nonzeros);
   }

   template <typename ElementType>
   void PrimalRegularization<ElementType>::initialize_statistics(Statistics& statistics, const Options& options) {
      statistics.add_column("regulariz", Statistics::double_width - 4, options.get_int("statistics_regularization_column_order"));
   }

   // Nocedal and Wright, p51
   template <typename ElementType>
   void PrimalRegularization<ElementType>::regularize_hessian(Statistics& statistics, SymmetricMatrix<size_t, ElementType>& hessian,
         const Inertia& expected_inertia) {
      DEBUG << "Current Hessian:\n" << hessian << '\n';
      const size_t hessian_dimension = hessian.dimension();
      const double smallest_diagonal_entry = hessian.smallest_diagonal_entry(hessian_dimension);
      DEBUG << "The minimal diagonal entry of the matrix is " << smallest_diagonal_entry << '\n';

      double regularization_factor = (smallest_diagonal_entry > 0.) ? 0. : this->regularization_initial_value - smallest_diagonal_entry;
      bool good_inertia = false;
      bool symbolic_analysis_performed = false;
      while (!good_inertia) {
         DEBUG << "Testing factorization with regularization factor " << regularization_factor << '\n';
         if (0. < regularization_factor) {
            hessian.set_regularization([=](size_t variable_index) {
               return (variable_index < hessian_dimension) ? regularization_factor : 0.;
            });
         }
         DEBUG << "Current Hessian:\n" << hessian << '\n';

         // perform the symbolic analysis only once
         if (!symbolic_analysis_performed) {
            this->linear_solver->do_symbolic_analysis(hessian);
            symbolic_analysis_performed = true;
         }
         this->linear_solver->do_numerical_factorization(hessian);
         const Inertia estimated_inertia = this->linear_solver->get_inertia();
         DEBUG << "Expected inertia: " << expected_inertia << '\n';
         DEBUG << "Estimated inertia: " << estimated_inertia << '\n';
         if (estimated_inertia == expected_inertia) {
            good_inertia = true;
            DEBUG << "Factorization was a success\n";
         }
         else {
            regularization_factor = (regularization_factor == 0.) ? this->regularization_initial_value : this->regularization_increase_factor * regularization_factor;
            if (regularization_factor > this->regularization_failure_threshold) {
               throw UnstableRegularization();
            }
         }
      }
      statistics.set("regulariz", regularization_factor);
   }

   template <typename ElementType>
   void PrimalRegularization<ElementType>::regularize_augmented_matrix(Statistics& /*statistics*/, SymmetricMatrix <size_t, ElementType>& /*augmented_matrix*/,
         const Collection<size_t>& /*primal_variables*/, const Collection<size_t>& /*dual_variables*/, ElementType /*dual_regularization_parameter*/,
         const Inertia& /*expected_inertia*/) {
      throw std::runtime_error("PrimalRegularization<ElementType>::regularize_augmented_matrix not implemented yet");
   }

} // namespace

#endif // UNO_PRIMALREGULARIZATION_H