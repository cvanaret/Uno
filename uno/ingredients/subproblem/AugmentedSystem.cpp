// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "AugmentedSystem.hpp"
#include "linear_algebra/SymmetricMatrixFactory.hpp"

AugmentedSystem::AugmentedSystem(const std::string& sparse_format, size_t max_dimension, size_t max_number_non_zeros, bool use_regularization,
      const Options& options):
   matrix(SymmetricMatrixFactory::create(sparse_format, max_dimension, max_number_non_zeros, use_regularization)),
   rhs(max_dimension),
   solution(max_dimension),
   regularization_failure_threshold(stod(options.at("regularization_failure_threshold"))),
   regularization_first_block_initial_factor(stod(options.at("regularization_first_block_initial_factor"))),
   regularization_second_block_fraction(stod(options.at("regularization_second_block_fraction"))),
   regularization_first_block_lb(stod(options.at("regularization_first_block_lb"))),
   regularization_first_block_decrease_factor(stod(options.at("regularization_first_block_decrease_factor"))),
   regularization_first_block_fast_increase_factor(stod(options.at("regularization_first_block_fast_increase_factor"))),
   regularization_first_block_slow_increase_factor(stod(options.at("regularization_first_block_slow_increase_factor"))) {
}

void AugmentedSystem::factorize_matrix(const Model& model, LinearSolver& linear_solver) {
   // compute the symbolic factorization only when:
   // the problem has a non-constant augmented system (ie is not an LP or a QP) or it is the first factorization
   if (this->number_factorizations == 0 || model.problem_type == NONLINEAR || !model.fixed_hessian_sparsity) {
      linear_solver.do_symbolic_factorization(*this->matrix);
   }
   linear_solver.do_numerical_factorization(*this->matrix);
   this->number_factorizations++;
}

void AugmentedSystem::regularize_matrix(const Model& model, LinearSolver& linear_solver, size_t size_first_block, size_t size_second_block,
      double constraint_regularization_parameter) {
   DEBUG << "Original matrix\n" << *this->matrix << '\n';
   double regularization_first_block = 0.;
   double regularization_second_block = 0.;
   DEBUG << "Testing factorization with regularization factors (" << regularization_first_block << ", " << regularization_second_block << ")\n";

   if (!linear_solver.matrix_is_singular() && linear_solver.number_negative_eigenvalues() == size_second_block) {
      DEBUG << "Inertia is good\n";
      return;
   }
   DEBUG << "Inertia is not good\n";
   // set the constraint regularization coefficient
   if (linear_solver.matrix_is_singular()) {
      DEBUG << "Matrix is singular\n";
      regularization_second_block = this->regularization_second_block_fraction * constraint_regularization_parameter;
   }
   else {
      regularization_second_block = 0.;
   }
   // set the Hessian regularization coefficient
   if (this->previous_regularization_first_block == 0.) {
      regularization_first_block = this->regularization_first_block_initial_factor;
   }
   else {
      regularization_first_block = std::max(this->regularization_first_block_lb,
            this->previous_regularization_first_block / this->regularization_first_block_decrease_factor);
   }

   // regularize the augmented matrix
   this->matrix->set_regularization([=](size_t i) {
      return (i < size_first_block) ? regularization_first_block : -regularization_second_block;
   });

   bool good_inertia = false;
   while (!good_inertia) {
      DEBUG << "Testing factorization with regularization factors (" << regularization_first_block << ", " << regularization_second_block << ")\n";
      DEBUG << *this->matrix << '\n';
      this->factorize_matrix(model, linear_solver);

      if (!linear_solver.matrix_is_singular() && linear_solver.number_negative_eigenvalues() == size_second_block) {
         good_inertia = true;
         DEBUG << "Factorization was a success\n\n";
         this->previous_regularization_first_block = regularization_first_block;
      }
      else {
         auto[number_pos_eigenvalues, number_neg_eigenvalues, number_zero_eigenvalues] = linear_solver.get_inertia();
         DEBUG << "Expected inertia (" << size_first_block << ", " << size_second_block << ", 0), ";
         DEBUG << "got (" << number_pos_eigenvalues << ", " << number_neg_eigenvalues << ", " << number_zero_eigenvalues << ")\n";
         if (this->previous_regularization_first_block == 0.) {
            regularization_first_block *= this->regularization_first_block_fast_increase_factor;
         }
         else {
            regularization_first_block *= this->regularization_first_block_slow_increase_factor;
         }

         if (regularization_first_block <= this->regularization_failure_threshold) {
            // regularize the augmented matrix
            this->matrix->set_regularization([=](size_t i) {
               return (i < size_first_block) ? regularization_first_block : -regularization_second_block;
            });
         }
         else {
            throw UnstableRegularization();
         }
      }
   }
}

void AugmentedSystem::solve(LinearSolver& linear_solver) {
   linear_solver.solve(*this->matrix, this->rhs, this->solution);
}