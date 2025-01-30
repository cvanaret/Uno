// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRIMALDUALCONVEXIFICATIONSTRATEGY_H
#define UNO_PRIMALDUALCONVEXIFICATIONSTRATEGY_H

#include <memory>
#include "ingredients/hessian_models/UnstableRegularization.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   // forward declarations
   template <typename IndexType, typename ElementType>
   class SymmetricMatrix;
   class WarmstartInformation;

   template <typename ElementType>
   class PrimalDualConvexificationStrategy {
   public:
      explicit PrimalDualConvexificationStrategy(const Options& options);
      void regularize_matrix(Statistics& statistics, DirectSymmetricIndefiniteLinearSolver<size_t, ElementType>& linear_solver,
            SymmetricMatrix<size_t, ElementType>& matrix, size_t size_primal_block, size_t size_dual_block, ElementType dual_regularization_parameter);

   protected:
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
   PrimalDualConvexificationStrategy<ElementType>::PrimalDualConvexificationStrategy(const Options& options):
         regularization_failure_threshold(ElementType(options.get_double("regularization_failure_threshold"))),
         primal_regularization_initial_factor(ElementType(options.get_double("primal_regularization_initial_factor"))),
         dual_regularization_fraction(ElementType(options.get_double("dual_regularization_fraction"))),
         primal_regularization_lb(ElementType(options.get_double("primal_regularization_lb"))),
         primal_regularization_decrease_factor(ElementType(options.get_double("primal_regularization_decrease_factor"))),
         primal_regularization_fast_increase_factor(ElementType(options.get_double("primal_regularization_fast_increase_factor"))),
         primal_regularization_slow_increase_factor(ElementType(options.get_double("primal_regularization_slow_increase_factor"))),
         threshold_unsuccessful_attempts(options.get_unsigned_int("threshold_unsuccessful_attempts")) {
   }

   // the matrix has been factorized prior to calling this function
   template <typename ElementType>
   void PrimalDualConvexificationStrategy<ElementType>::regularize_matrix(Statistics& statistics,
         DirectSymmetricIndefiniteLinearSolver<size_t, ElementType>& linear_solver, SymmetricMatrix<size_t, ElementType>& matrix,
         size_t size_primal_block, size_t size_dual_block, ElementType dual_regularization_parameter) {
      DEBUG2 << "Original matrix\n" << matrix << '\n';
      this->primal_regularization = ElementType(0.);
      this->dual_regularization = ElementType(0.);
      DEBUG << "Testing factorization with regularization factors (0, 0)\n";
      size_t number_attempts = 1;
      DEBUG << "Number of attempts: " << number_attempts << "\n\n";

      auto [number_pos_eigenvalues, number_neg_eigenvalues, number_zero_eigenvalues] = linear_solver.get_inertia();
      DEBUG << "Expected inertia  (" << size_primal_block << ", " << size_dual_block << ", 0)\n";
      DEBUG << "Estimated inertia (" << number_pos_eigenvalues << ", " << number_neg_eigenvalues << ", " << number_zero_eigenvalues << ")\n";

      if (number_pos_eigenvalues == size_primal_block && number_neg_eigenvalues == size_dual_block && number_zero_eigenvalues == 0) {
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
      matrix.set_regularization([=](size_t row_index) {
         return (row_index < size_primal_block) ? this->primal_regularization : -this->dual_regularization;
      });

      bool good_inertia = false;
      while (!good_inertia) {
         DEBUG << "Testing factorization with regularization factors (" << this->primal_regularization << ", " << this->dual_regularization << ")\n";
         DEBUG2 << matrix << '\n';
         DEBUG << "Performing numerical factorization of the indefinite system\n";
         linear_solver.do_numerical_factorization(matrix);
         number_attempts++;
         DEBUG << "Number of attempts: " << number_attempts << "\n";

         std::tie(number_pos_eigenvalues, number_neg_eigenvalues, number_zero_eigenvalues) = linear_solver.get_inertia();
         DEBUG << "Expected inertia  (" << size_primal_block << ", " << size_dual_block << ", 0)\n";
         DEBUG << "Estimated inertia (" << number_pos_eigenvalues << ", " << number_neg_eigenvalues << ", " << number_zero_eigenvalues << ")\n";

         if (number_pos_eigenvalues == size_primal_block && number_neg_eigenvalues == size_dual_block && number_zero_eigenvalues == 0) {
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
               matrix.set_regularization([=](size_t row_index) {
                  return (row_index < size_primal_block) ? this->primal_regularization : -this->dual_regularization;
               });
            }
            else {
               throw UnstableRegularization();
            }
         }
      }
      // check the inertia
      std::tie(number_pos_eigenvalues, number_neg_eigenvalues, number_zero_eigenvalues) = linear_solver.get_inertia();
      assert(number_pos_eigenvalues == size_primal_block && number_neg_eigenvalues == size_dual_block && number_zero_eigenvalues == 0);
      statistics.set("regulariz", this->primal_regularization);
   }
} // namespace

#endif // UNO_PRIMALDUALCONVEXIFICATIONSTRATEGY_H
