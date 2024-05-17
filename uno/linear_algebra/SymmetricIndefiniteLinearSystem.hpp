// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICINDEFINITELINEARSYSTEM_H
#define UNO_SYMMETRICINDEFINITELINEARSYSTEM_H

#include <memory>
#include "SymmetricMatrix.hpp"
#include "SymmetricMatrixFactory.hpp"
#include "RectangularMatrix.hpp"
#include "model/Model.hpp"
#include "solvers/linear/direct/DirectIndefiniteLinearSolver.hpp"
#include "tools/Options.hpp"
#include "tools/Statistics.hpp"

struct UnstableRegularization : public std::exception {

   [[nodiscard]] const char* what() const noexcept override {
      return "The inertia correction got unstable (delta_w > threshold)";
   }
};

template <typename IndexType, typename ElementType>
class SymmetricIndefiniteLinearSystem {
public:
   std::unique_ptr<SymmetricMatrix<IndexType, ElementType>> matrix;
   std::vector<ElementType> rhs{};
   std::vector<ElementType> solution{};

   SymmetricIndefiniteLinearSystem(const std::string& sparse_format, size_t dimension, size_t number_non_zeros, bool use_regularization,
         const Options& options);
   void assemble_matrix(const SymmetricMatrix<IndexType, ElementType>& hessian, const RectangularMatrix<ElementType>& constraint_jacobian,
         size_t number_variables, size_t number_constraints);
   void factorize_matrix(DirectIndefiniteLinearSolver<IndexType, ElementType>& linear_solver, bool fixed_sparsity_pattern);
   void regularize_matrix(Statistics& statistics, DirectIndefiniteLinearSolver<IndexType, ElementType>& linear_solver, bool fixed_sparsity_pattern,
         size_t size_primal_block, size_t size_dual_block, ElementType dual_regularization_parameter);
   void solve(DirectIndefiniteLinearSolver<IndexType, ElementType>& linear_solver);

protected:
   size_t number_factorizations{0};
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

template <typename IndexType, typename ElementType>
SymmetricIndefiniteLinearSystem<IndexType, ElementType>::SymmetricIndefiniteLinearSystem(const std::string& sparse_format, size_t dimension,
      size_t number_non_zeros, bool use_regularization, const Options& options):
      matrix(SymmetricMatrixFactory<IndexType, ElementType>::create(sparse_format, dimension, number_non_zeros, use_regularization)),
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

template <typename IndexType, typename ElementType>
void SymmetricIndefiniteLinearSystem<IndexType, ElementType>::assemble_matrix(const SymmetricMatrix<IndexType, ElementType>& hessian,
      const RectangularMatrix<ElementType>& constraint_jacobian, size_t number_variables, size_t number_constraints) {
   this->matrix->dimension = number_variables + number_constraints;
   this->matrix->reset();
   // copy the Lagrangian Hessian in the top left block
   //size_t current_column = 0;
   hessian.for_each([&](IndexType row_index, IndexType column_index, ElementType entry) {
      // finalize all empty columns
      /*for (size_t column: Range(current_column, column_index)) {
         this->matrix->finalize_column(column);
         current_column++;
      }*/
      this->matrix->insert(entry, row_index, column_index);
   });

   // Jacobian of general constraints
   for (size_t column_index: Range(number_constraints)) {
      constraint_jacobian[column_index].for_each([&](IndexType row_index, ElementType derivative) {
         this->matrix->insert(derivative, row_index, number_variables + column_index);
      });
      this->matrix->finalize_column(column_index);
   }
}

template <typename IndexType, typename ElementType>
void SymmetricIndefiniteLinearSystem<IndexType, ElementType>::factorize_matrix(DirectIndefiniteLinearSolver<IndexType, ElementType>& linear_solver, bool fixed_sparsity_pattern) {
   // compute the symbolic factorization only when:
   // the problem has a non-constant augmented system (ie is not an LP or a QP) or it is the first factorization
   if (true || this->number_factorizations == 0 || not fixed_sparsity_pattern) {
      linear_solver.do_symbolic_factorization(*this->matrix);
   }
   linear_solver.do_numerical_factorization(*this->matrix);
   this->number_factorizations++;
}

template <typename IndexType, typename ElementType>
void SymmetricIndefiniteLinearSystem<IndexType, ElementType>::regularize_matrix(Statistics& statistics, DirectIndefiniteLinearSolver<IndexType, ElementType>& linear_solver,
      bool fixed_sparsity_pattern, size_t size_primal_block, size_t size_dual_block, ElementType dual_regularization_parameter) {
   DEBUG2 << "Original matrix\n" << *this->matrix << '\n';
   this->primal_regularization = ElementType(0.);
   this->dual_regularization = ElementType(0.);
   size_t number_attempts = 1;
   DEBUG << "Testing factorization with regularization factors (" << this->primal_regularization << ", " << this->dual_regularization << ")\n";

   if (not linear_solver.matrix_is_singular() && linear_solver.number_negative_eigenvalues() == size_dual_block) {
      DEBUG << "Inertia is good\n";
      statistics.set("regularization", this->primal_regularization);
      return;
   }
   auto[number_pos_eigenvalues, number_neg_eigenvalues, number_zero_eigenvalues] = linear_solver.get_inertia();
   DEBUG << "Expected inertia (" << size_primal_block << ", " << size_dual_block << ", 0), ";
   DEBUG << "got (" << number_pos_eigenvalues << ", " << number_neg_eigenvalues << ", " << number_zero_eigenvalues << ")\n";
   DEBUG << "Number of attempts: " << number_attempts << "\n\n";

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
      this->primal_regularization = std::max(this->primal_regularization_lb, this->previous_primal_regularization / this->primal_regularization_decrease_factor);
   }

   // regularize the augmented matrix
   this->matrix->set_regularization([=](size_t row_index) {
      return (row_index < size_primal_block) ? this->primal_regularization : -this->dual_regularization;
   });

   bool good_inertia = false;
   while (not good_inertia) {
      DEBUG << "Testing factorization with regularization factors (" << this->primal_regularization << ", " << this->dual_regularization << ")\n";
      DEBUG2 << *this->matrix << '\n';
      this->factorize_matrix(linear_solver, fixed_sparsity_pattern);
      number_attempts++;

      if (not linear_solver.matrix_is_singular() && linear_solver.number_negative_eigenvalues() == size_dual_block) {
         good_inertia = true;
         DEBUG << "Factorization was a success\n\n";
         this->previous_primal_regularization = this->primal_regularization;
      }
      else {
         auto[number_pos_eigenvalues, number_neg_eigenvalues, number_zero_eigenvalues] = linear_solver.get_inertia();
         DEBUG << "Expected inertia (" << size_primal_block << ", " << size_dual_block << ", 0), ";
         DEBUG << "got (" << number_pos_eigenvalues << ", " << number_neg_eigenvalues << ", " << number_zero_eigenvalues << ")\n";
         DEBUG << "Number of attempts: " << number_attempts << "\n";
         if (this->previous_primal_regularization == 0. || this->threshold_unsuccessful_attempts < number_attempts) {
            this->primal_regularization *= this->primal_regularization_fast_increase_factor;
         }
         else {
            this->primal_regularization *= this->primal_regularization_slow_increase_factor;
         }

         if (this->primal_regularization <= this->regularization_failure_threshold) {
            // regularize the augmented matrix
            this->matrix->set_regularization([=](size_t row_index) {
               return (row_index < size_primal_block) ? this->primal_regularization : -this->dual_regularization;
            });
         }
         else {
            throw UnstableRegularization();
         }
      }
   }
   statistics.set("regularization", this->primal_regularization);
}

template <typename IndexType, typename ElementType>
void SymmetricIndefiniteLinearSystem<IndexType, ElementType>::solve(DirectIndefiniteLinearSolver<IndexType, ElementType>& linear_solver) {
   const bool from_scratch = (this->number_factorizations == 0);
   linear_solver.solve_indefinite_system(*this->matrix, this->rhs, this->solution, from_scratch);
}

/*
template <typename IndexType, typename ElementType>
ElementType SymmetricIndefiniteLinearSystem<IndexType, ElementType>::get_primal_regularization() const {
   return this->primal_regularization;
}
*/

#endif // UNO_SYMMETRICINDEFINITELINEARSYSTEM_H
