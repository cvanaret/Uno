// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRIMALREGULARIZATION_H
#define UNO_PRIMALREGULARIZATION_H

#include "UnstableRegularization.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   template <typename ElementType>
   class PrimalRegularization: public RegularizationStrategy<ElementType> {
   public:
      explicit PrimalRegularization(const Options& options);

       void initialize_statistics(Statistics& statistics, const Options& options) override;

      void regularize_matrix(Statistics& statistics, DirectSymmetricIndefiniteLinearSolver<size_t, ElementType>& linear_solver,
            SymmetricMatrix<size_t, ElementType>& matrix, size_t size_primal_block, size_t size_dual_block,
            ElementType dual_regularization_parameter) override;

   protected:
      const double regularization_initial_value{};
      const double regularization_increase_factor{};
      const double regularization_failure_threshold{};
   };

   template <typename ElementType>
   PrimalRegularization<ElementType>::PrimalRegularization(const Options& options):
         RegularizationStrategy<ElementType>(),
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
   void PrimalRegularization<ElementType>::regularize_matrix(Statistics& statistics,
         DirectSymmetricIndefiniteLinearSolver <size_t, ElementType>& linear_solver, SymmetricMatrix <size_t, ElementType>& matrix,
         size_t size_primal_block, size_t /*size_dual_block*/, ElementType /*dual_regularization_parameter*/) {
      DEBUG << "Current matrix:\n" << matrix << '\n';
      const double smallest_diagonal_entry = matrix.smallest_diagonal_entry(size_primal_block); // TODO check that
      DEBUG << "The minimal diagonal entry of the matrix is " << smallest_diagonal_entry << '\n';

      double regularization_factor = (smallest_diagonal_entry > 0.) ? 0. : this->regularization_initial_value - smallest_diagonal_entry;
      bool good_inertia = false;
      bool symbolic_analysis_performed = false;
      while (!good_inertia) {
         DEBUG << "Testing factorization with regularization factor " << regularization_factor << '\n';
         if (0. < regularization_factor) {
            matrix.set_regularization([=](size_t variable_index) {
               return (variable_index < size_primal_block) ? regularization_factor : 0.;
            });
         }
         DEBUG << "Current matrix:\n" << matrix << '\n';

         // perform the symbolic analysis only once
         if (!symbolic_analysis_performed) {
            linear_solver.do_symbolic_analysis(matrix);
            symbolic_analysis_performed = true;
         }
         linear_solver.do_numerical_factorization(matrix);
         if (linear_solver.rank() == size_primal_block && linear_solver.number_negative_eigenvalues() == 0) {
            good_inertia = true;
            DEBUG << "Factorization was a success\n";
         }
         else {
            DEBUG << "rank: " << linear_solver.rank() << ", negative eigenvalues: " << linear_solver.number_negative_eigenvalues() << '\n';
            regularization_factor = (regularization_factor == 0.) ? this->regularization_initial_value : this->regularization_increase_factor * regularization_factor;
            if (regularization_factor > this->regularization_failure_threshold) {
               throw UnstableRegularization();
            }
         }
      }
      statistics.set("regulariz", regularization_factor);
   }

} // namespace

#endif // UNO_PRIMALREGULARIZATION_H