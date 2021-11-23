#include "AugmentedSystem.hpp"
#include "linear_algebra/SymmetricMatrixFactory.hpp"

AugmentedSystem::AugmentedSystem(const std::string& sparse_format, size_t max_dimension, size_t max_number_non_zeros):
   matrix(SymmetricMatrixFactory::create(sparse_format, max_dimension, max_number_non_zeros)),
   rhs(max_dimension), solution(max_dimension) {
}

void AugmentedSystem::solve(LinearSolver& linear_solver, size_t dimension) {
   linear_solver.solve(dimension, *this->matrix, this->rhs, this->solution);
}

void AugmentedSystem::factorize_matrix(const Problem& problem, LinearSolver& linear_solver, size_t dimension) {
   // compute the symbolic factorization only when:
   // the problem has a non-constant augmented system (ie is not an LP or a QP) or it is the first factorization
   if (this->number_factorizations == 0 || !problem.fixed_hessian_sparsity || problem.problem_type == NONLINEAR) {
      linear_solver.do_symbolic_factorization(dimension, *this->matrix);
   }
   linear_solver.do_numerical_factorization(dimension, *this->matrix);
   this->number_factorizations++;
}

void AugmentedSystem::regularize_matrix(const Problem& problem, LinearSolver& linear_solver, size_t size_first_block, size_t size_second_block,
      double constraint_regularization_parameter) {
   DEBUG << "Original matrix\n" << *this->matrix << "\n";
   this->regularization_first_block = 0.;
   this->regularization_second_block = 0.;
   DEBUG << "Testing factorization with regularization factor " << this->regularization_first_block << "\n";

   bool good_inertia = false;
   if (!linear_solver.matrix_is_singular() && linear_solver.number_negative_eigenvalues() == size_second_block) {
      DEBUG << "Inertia is good\n";
      good_inertia = true;
   }
   else {
      DEBUG << "Inertia is not good\n";
      // constraint regularization
      if (linear_solver.matrix_is_singular()) {
         DEBUG << "Matrix is singular\n";
         this->regularization_second_block = 1e-8 * constraint_regularization_parameter;
      }
      else {
         this->regularization_second_block = 0.;
      }
      // Hessian regularization
      if (this->previous_regularization_first_block == 0.) {
         this->regularization_first_block = 1e-4;
      }
      else {
         this->regularization_first_block = std::max(1e-20, this->previous_regularization_first_block / 3.);
      }
   }

   size_t current_matrix_size = this->matrix->number_nonzeros;
   if (!good_inertia) {
      for (size_t i = 0; i < size_first_block; i++) {
         this->matrix->insert(this->regularization_first_block, i, i);
      }
      for (size_t j = size_first_block; j < size_first_block + size_second_block; j++) {
         this->matrix->insert(-this->regularization_second_block, j, j);
      }
   }

   while (!good_inertia) {
      DEBUG << "Testing factorization with regularization factor " << this->regularization_first_block << "\n";
      DEBUG << *this->matrix << "\n";
      this->factorize_matrix(problem, linear_solver, size_first_block + size_second_block);

      if (!linear_solver.matrix_is_singular() && linear_solver.number_negative_eigenvalues() == size_second_block) {
         good_inertia = true;
         DEBUG << "Factorization was a success\n\n";
         this->previous_regularization_first_block = this->regularization_first_block;
      }
      else {
         if (this->previous_regularization_first_block == 0.) {
            this->regularization_first_block *= 100.;
         }
         else {
            this->regularization_first_block *= 8.;
         }

         if (this->regularization_first_block <= this->regularization_failure_threshold) {
            for (size_t i = 0; i < size_first_block; i++) {
               this->matrix->entries[current_matrix_size + i] = this->regularization_first_block;
            }
            for (size_t j = size_first_block; j < size_first_block + size_second_block; j++) {
               this->matrix->entries[current_matrix_size + j] = -this->regularization_second_block;
            }
         }
         else {
            throw UnstableInertiaCorrection();
         }
      }
   }
}