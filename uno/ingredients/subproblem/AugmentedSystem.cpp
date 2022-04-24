#include "AugmentedSystem.hpp"
#include "linear_algebra/SymmetricMatrixFactory.hpp"

AugmentedSystem::AugmentedSystem(const std::string& sparse_format, size_t max_dimension, size_t max_number_non_zeros,
      double regularization_failure_threshold):
   matrix(SymmetricMatrixFactory::create(sparse_format, max_dimension, max_number_non_zeros)),
   rhs(max_dimension),
   solution(max_dimension),
   regularization_failure_threshold(regularization_failure_threshold) {
}

void AugmentedSystem::solve(LinearSolver& linear_solver) {
   linear_solver.solve(*this->matrix, this->rhs, this->solution);
}

void AugmentedSystem::factorize_matrix(const ReformulatedProblem& /*problem*/, LinearSolver& linear_solver) {
   // compute the symbolic factorization only when:
   // the problem has a non-constant augmented system (ie is not an LP or a QP) or it is the first factorization
   if (this->number_factorizations == 0 || true) { // || !model.fixed_hessian_sparsity || model.problem_type == NONLINEAR) {
      // TODO WARNING << "AugmentedSystem::factorize_matrix: handle fixed_hessian_sparsity and problem_type\n";
      linear_solver.do_symbolic_factorization(*this->matrix);
   }
   linear_solver.do_numerical_factorization(*this->matrix);
   this->number_factorizations++;
}

void AugmentedSystem::regularize_matrix(const ReformulatedProblem& problem, LinearSolver& linear_solver, size_t size_first_block, size_t size_second_block,
      double constraint_regularization_parameter) {
   DEBUG << "Original matrix\n" << *this->matrix << '\n';
   double regularization_first_block = 0.;
   double regularization_second_block = 0.;
   DEBUG << "Testing factorization with regularization factor " << regularization_first_block << '\n';

   if (!linear_solver.matrix_is_singular() && linear_solver.number_negative_eigenvalues() == size_second_block) {
      DEBUG << "Inertia is good\n";
      return;
   }
   DEBUG << "Inertia is not good\n";
   // set the constraint regularization coefficient
   if (linear_solver.matrix_is_singular()) {
      DEBUG << "Matrix is singular\n";
      regularization_second_block = 1e-8 * constraint_regularization_parameter;
   }
   else {
      regularization_second_block = 0.;
   }
   // set the Hessian regularization coefficient
   if (this->previous_regularization_first_block == 0.) {
      regularization_first_block = 1e-4;
   }
   else {
      regularization_first_block = std::max(1e-20, this->previous_regularization_first_block / 3.);
   }
   size_t current_matrix_size = this->matrix->number_nonzeros;

   // regularize
   for (size_t i = 0; i < size_first_block; i++) {
      this->matrix->insert(regularization_first_block, i, i);
   }
   for (size_t j = size_first_block; j < size_first_block + size_second_block; j++) {
      this->matrix->insert(-regularization_second_block, j, j);
   }

   bool good_inertia = false;
   while (!good_inertia) {
      DEBUG << "Testing factorization with regularization factor " << regularization_first_block << '\n';
      DEBUG << *this->matrix << '\n';
      this->factorize_matrix(problem, linear_solver);

      if (!linear_solver.matrix_is_singular() && linear_solver.number_negative_eigenvalues() == size_second_block) {
         good_inertia = true;
         DEBUG << "Factorization was a success\n\n";
         this->previous_regularization_first_block = regularization_first_block;
      }
      else {
         if (this->previous_regularization_first_block == 0.) {
            regularization_first_block *= 100.;
         }
         else {
            regularization_first_block *= 8.;
         }

         if (regularization_first_block <= this->regularization_failure_threshold) {
            for (size_t i = 0; i < size_first_block; i++) {
               this->matrix->entries[current_matrix_size + i] = regularization_first_block;
            }
            for (size_t j = size_first_block; j < size_first_block + size_second_block; j++) {
               this->matrix->entries[current_matrix_size + j] = -regularization_second_block;
            }
         }
         else {
            throw UnstableRegularization();
         }
      }
   }
}
