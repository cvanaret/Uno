// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICINDEFINITELINEARSYSTEM_H
#define UNO_SYMMETRICINDEFINITELINEARSYSTEM_H

#include "SymmetricMatrix.hpp"
#include "SparseStorageFactory.hpp"
#include "RectangularMatrix.hpp"
#include "../ingredients/regularization_strategies/UnstableRegularization.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   template <typename ElementType>
   class SymmetricIndefiniteLinearSystem {
   public:
      SymmetricMatrix<size_t, ElementType> matrix{};
      Vector<ElementType> rhs{};
      Vector<ElementType> solution{};

      SymmetricIndefiniteLinearSystem(const std::string& sparse_format, size_t dimension, size_t number_non_zeros, bool use_regularization,
            const Options& options);
      explicit SymmetricIndefiniteLinearSystem(const Options& options);
      //SymmetricIndefiniteLinearSystem& operator=(const SymmetricIndefiniteLinearSystem& other) = default;
      SymmetricIndefiniteLinearSystem& operator=(SymmetricIndefiniteLinearSystem&& other) = default;

      void assemble_matrix(const SymmetricMatrix<size_t, double>& hessian, const RectangularMatrix<double>& constraint_jacobian,
            size_t number_variables, size_t number_constraints);
      void factorize_matrix(DirectSymmetricIndefiniteLinearSolver<size_t, ElementType>& linear_solver, WarmstartInformation& warmstart_information);
      void regularize_matrix(Statistics& statistics, DirectSymmetricIndefiniteLinearSolver<size_t, ElementType>& linear_solver,
            size_t size_primal_block, size_t size_dual_block, ElementType dual_regularization_parameter, WarmstartInformation& warmstart_information);
      void solve(DirectSymmetricIndefiniteLinearSolver<size_t, ElementType>& linear_solver);
      // [[nodiscard]] T get_primal_regularization() const;

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
   SymmetricIndefiniteLinearSystem<ElementType>::SymmetricIndefiniteLinearSystem(const std::string& sparse_format, size_t dimension,
         size_t number_non_zeros, bool use_regularization, const Options& options):
         matrix(dimension, number_non_zeros, use_regularization, sparse_format),
         rhs(dimension),
         solution(dimension),
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
   SymmetricIndefiniteLinearSystem<ElementType>::SymmetricIndefiniteLinearSystem(const Options& options):
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
   void SymmetricIndefiniteLinearSystem<ElementType>::assemble_matrix(const SymmetricMatrix<size_t, double>& hessian,
         const RectangularMatrix<double>& constraint_jacobian, size_t number_variables, size_t number_constraints) {
      this->matrix.set_dimension(number_variables + number_constraints);
      this->matrix.reset();
      // copy the Lagrangian Hessian in the top left block
      //size_t current_column = 0;
      for (const auto [row_index, column_index, element]: hessian) {
         // finalize all empty columns
         /*for (size_t column: Range(current_column, column_index)) {
            this->matrix.finalize_column(column);
            current_column++;
         }*/
         this->matrix.insert(element, row_index, column_index);
      }

      // Jacobian of general constraints
      for (size_t column_index: Range(number_constraints)) {
         for (const auto [row_index, derivative]: constraint_jacobian[column_index]) {
            this->matrix.insert(derivative, row_index, number_variables + column_index);
         }
         this->matrix.finalize_column(column_index);
      }
   }

   template <typename ElementType>
   void SymmetricIndefiniteLinearSystem<ElementType>::factorize_matrix(DirectSymmetricIndefiniteLinearSolver<size_t, ElementType>& linear_solver,
         WarmstartInformation& warmstart_information) {
      if (warmstart_information.hessian_sparsity_changed || warmstart_information.jacobian_sparsity_changed) {
         DEBUG << "Performing symbolic analysis of the indefinite system\n";
         linear_solver.do_symbolic_analysis(this->matrix);
         warmstart_information.hessian_sparsity_changed = warmstart_information.jacobian_sparsity_changed = false;
      }
      DEBUG << "Performing numerical factorization of the indefinite system\n";
      linear_solver.do_numerical_factorization(this->matrix);
   }

   // the matrix has been factorized prior to calling this function
   template <typename ElementType>
   void SymmetricIndefiniteLinearSystem<ElementType>::regularize_matrix(Statistics& statistics,
         DirectSymmetricIndefiniteLinearSolver<size_t, ElementType>& linear_solver, size_t size_primal_block, size_t size_dual_block,
         ElementType dual_regularization_parameter, WarmstartInformation& warmstart_information) {
      DEBUG2 << "Original matrix\n" << this->matrix << '\n';
      this->primal_regularization = ElementType(0.);
      this->dual_regularization = ElementType(0.);
      size_t number_attempts = 1;
      DEBUG << "Number of attempts: " << number_attempts << "\n\n";

      const Inertia estimated_inertia = linear_solver.get_inertia();
      DEBUG << "Expected inertia  (" << size_primal_block << ", " << size_dual_block << ", 0)\n";
      DEBUG << "Estimated inertia " << estimated_inertia << '\n';

      if (estimated_inertia.positive == size_primal_block && estimated_inertia.negative == size_dual_block && estimated_inertia.zero == 0) {
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
      this->matrix.set_regularization([=](size_t row_index) {
         return (row_index < size_primal_block) ? this->primal_regularization : -this->dual_regularization;
      });

      bool good_inertia = false;
      while (!good_inertia) {
         DEBUG << "Testing factorization with regularization factors (" << this->primal_regularization << ", " << this->dual_regularization << ")\n";
         DEBUG2 << this->matrix << '\n';
         this->factorize_matrix(linear_solver, warmstart_information);
         number_attempts++;
         DEBUG << "Number of attempts: " << number_attempts << "\n";

         const Inertia new_estimated_inertia = linear_solver.get_inertia();
         DEBUG << "Expected inertia  (" << size_primal_block << ", " << size_dual_block << ", 0)\n";
         DEBUG << "Estimated inertia " << new_estimated_inertia << '\n';

         if (new_estimated_inertia.positive == size_primal_block && new_estimated_inertia.negative == size_dual_block && new_estimated_inertia.zero == 0) {
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
               this->matrix.set_regularization([=](size_t row_index) {
                  return (row_index < size_primal_block) ? this->primal_regularization : -this->dual_regularization;
               });
            }
            else {
               throw UnstableRegularization();
            }
         }
      }
      statistics.set("regulariz", this->primal_regularization);
   }

   template <typename ElementType>
   void SymmetricIndefiniteLinearSystem<ElementType>::solve(DirectSymmetricIndefiniteLinearSolver<size_t, ElementType>& linear_solver) {
      linear_solver.solve_indefinite_system(this->matrix, this->rhs, this->solution);
   }

   /*
   template <typename ElementType>
   ElementType SymmetricIndefiniteLinearSystem<ElementType>::get_primal_regularization() const {
      return this->primal_regularization;
   }
   */
} // namespace

#endif // UNO_SYMMETRICINDEFINITELINEARSYSTEM_H
