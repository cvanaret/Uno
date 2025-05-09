// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRIMALDUALREGULARIZATION_H
#define UNO_PRIMALDUALREGULARIZATION_H

#include <cassert>
#include "RegularizationStrategy.hpp"
#include "UnstableRegularization.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   // forward declaration
   template <typename IndexType, typename ElementType>
   class SymmetricMatrix;

   template <typename ElementType>
   class PrimalDualRegularization: public RegularizationStrategy<ElementType> {
   public:
      explicit PrimalDualRegularization(const Options& options);

      void initialize_memory(const OptimizationProblem& problem, const HessianModel& hessian_model) override;
      void initialize_statistics(Statistics& statistics, const Options& options) override;

      void regularize_hessian(Statistics& statistics, SymmetricMatrix<size_t, ElementType>& hessian, const Inertia& expected_inertia) override;
      void regularize_augmented_matrix(Statistics& statistics, SymmetricMatrix<size_t, ElementType>& augmented_matrix, const Collection<size_t>& primal_variables,
         const Collection<size_t>& dual_variables, ElementType dual_regularization_parameter, const Inertia& expected_inertia) override;

   protected:
      std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<size_t, double>> linear_solver;
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
         linear_solver(SymmetricIndefiniteLinearSolverFactory::create(options)),
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
   void PrimalDualRegularization<ElementType>::initialize_memory(const OptimizationProblem& problem, const HessianModel& hessian_model) {
      const size_t dimension = problem.number_variables + problem.number_constraints;
      const size_t number_nonzeros = problem.number_hessian_nonzeros(hessian_model) + problem.number_jacobian_nonzeros() +
         problem.number_variables + problem.number_constraints; // diagonal primal-dual regularization
      this->linear_solver->initialize_memory(dimension, number_nonzeros);
   }

   template <typename ElementType>
   void PrimalDualRegularization<ElementType>::initialize_statistics(Statistics& statistics, const Options& options) {
      statistics.add_column("regulariz", Statistics::double_width - 4, options.get_int("statistics_regularization_column_order"));
   }

   template <typename ElementType>
   void PrimalDualRegularization<ElementType>::regularize_hessian(Statistics& statistics, SymmetricMatrix<size_t, ElementType>& hessian,
         const Inertia& expected_inertia) {
      // to regularize the Hessian only, call the function for the augmented matrix with no dual part
      this->regularize_augmented_matrix(statistics, hessian, Range(hessian.dimension()), Range(0), ElementType(0), expected_inertia);
   }

   // the augmented matrix has been factorized prior to calling this function
   template <typename ElementType>
   void PrimalDualRegularization<ElementType>::regularize_augmented_matrix(Statistics& statistics, SymmetricMatrix<size_t, ElementType>& augmented_matrix,
         const Collection<size_t>& /*primal_variables*/, const Collection<size_t>& /*dual_variables*/, ElementType dual_regularization_parameter,
         const Inertia& expected_inertia) {
      this->primal_regularization = ElementType(0);
      this->dual_regularization = ElementType(0);
      DEBUG << "Testing factorization with regularization factors (0, 0)\n";
      size_t number_attempts = 1;
      DEBUG << "Number of attempts: " << number_attempts << "\n\n";

      const Inertia estimated_inertia = this->linear_solver->get_inertia();
      DEBUG << "Expected inertia  (" << expected_inertia.positive << ", " << expected_inertia.negative << ", " << expected_inertia.zero << ")\n";
      DEBUG << "Estimated inertia (" << estimated_inertia.positive << ", " << estimated_inertia.negative << ", " << estimated_inertia.zero << ")\n";

      if (estimated_inertia == expected_inertia) {
         DEBUG << "The inertia is correct\n";
         statistics.set("regulariz", this->primal_regularization);
         return;
      }

      // set the constraint regularization coefficient
      if (this->linear_solver->matrix_is_singular()) {
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
      augmented_matrix.set_regularization([=](size_t row_index) {
         return (row_index < expected_inertia.positive) ? this->primal_regularization : -this->dual_regularization; // TODO
      });

      bool good_inertia = false;
      while (!good_inertia) {
         DEBUG << "Testing factorization with regularization factors (" << this->primal_regularization << ", " << this->dual_regularization << ")\n";
         DEBUG2 << augmented_matrix << '\n';
         DEBUG << "Performing numerical factorization of the indefinite system\n";
         this->linear_solver->do_numerical_factorization(augmented_matrix);
         number_attempts++;
         DEBUG << "Number of attempts: " << number_attempts << "\n";

         const Inertia new_estimated_inertia = this->linear_solver->get_inertia();
         DEBUG << "Expected inertia  " << expected_inertia << ")\n";
         DEBUG << "Estimated inertia " << new_estimated_inertia << ")\n";

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
               augmented_matrix.set_regularization([=](size_t row_index) {
                  return (row_index < expected_inertia.positive) ? this->primal_regularization : -this->dual_regularization; // TODO
               });
            }
            else {
               throw UnstableRegularization();
            }
         }
      }
      // check the inertia
      const Inertia latest_estimated_inertia = this->linear_solver->get_inertia();
      assert(latest_estimated_inertia == expected_inertia);
      statistics.set("regulariz", this->primal_regularization);
   }
} // namespace

#endif // UNO_PRIMALDUALREGULARIZATION_H