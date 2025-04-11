// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include "LBFGSHessian.hpp"

#include <model/Model.hpp>

#include "linear_algebra/LAPACK.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/Range.hpp"

namespace uno {
   LBFGSHessian::LBFGSHessian(size_t number_variables, size_t memory_size):
         HessianModel(),
         number_variables(number_variables),
         memory_size(memory_size),
         current_memory_slot(memory_size-1),
         S_matrix(number_variables, memory_size),
         Y_matrix(number_variables, memory_size),
         L_matrix(memory_size, memory_size),
         D_matrix(memory_size),
         M_matrix(memory_size, memory_size) {
   }

   void LBFGSHessian::initialize_statistics(Statistics& /*statistics*/, const Options& /*options*/) const {
      // do nothing
   }

   void LBFGSHessian::notify_accepted_iterate(const Model& model, const Iterate& current_iterate, const Iterate& trial_iterate) {
      this->update_memory(model, current_iterate, trial_iterate);
      this->hessian_recomputation_required = true;
   }

   void LBFGSHessian::evaluate_hessian(Statistics& /*statistics*/, const Model& /*model*/, const Vector<double>& /*primal_variables*/,
         double /*objective_multiplier*/, const Vector<double>& /*constraint_multipliers*/, SymmetricMatrix<size_t, double>& /*hessian*/) {
      if (this->hessian_recomputation_required) {
         this->recompute_hessian_representation();
         this->hessian_recomputation_required = false;
      }

      // TODO
      throw std::runtime_error("LBFGSHessian::evaluate_hessian not implemented");
   }

   void LBFGSHessian::compute_hessian_vector_product(const Model& /*model*/, const Vector<double>& vector, double /*objective_multiplier*/,
         const Vector<double>& /*constraint_multipliers*/, Vector<double>& result) {
      if (this->hessian_recomputation_required) {
         this->recompute_hessian_representation();
         this->hessian_recomputation_required = false;
      }

      // for the moment, pretend we have an identity Hessian TODO
      for (size_t variable_index: Range(vector.size())) {
         result[variable_index] = vector[variable_index];
      }
   }

   void LBFGSHessian::update_memory(const Model& model, const Iterate& current_iterate, const Iterate& trial_iterate) {
      // increment the slot: if we exceed the size of the memory, we start over and replace the older point in memory
      this->current_memory_slot = (this->current_memory_slot + 1) % this->memory_size;
      std::cout << "\n*** Adding vector to L-BFGS memory at slot " << this->current_memory_slot << '\n';
      // TODO figure out if we're extending or replacing in memory

      // fill the S matrix
      for (size_t variable_index: Range(this->number_variables)) {
         this->S_matrix.entry(variable_index, this->current_memory_slot) = trial_iterate.primals[variable_index] -
            current_iterate.primals[variable_index];
      }
      std::cout << "S:\n" << this->S_matrix;

      // fill the Y matrix
      // TODO
      std::cout << "Current Lag gradient: " << current_iterate.residuals.lagrangian_gradient.assemble(1.);
      std::cout << "Trial Lag gradient: " << trial_iterate.residuals.lagrangian_gradient.assemble(1.);

      this->used_memory_size = std::min(this->used_memory_size + 1, this->memory_size);
      std::cout << "There are now " << this->used_memory_size << " iterates in memory (capacity " <<
         this->memory_size << ")\n";
   }

   void LBFGSHessian::recompute_hessian_representation() {
      std::cout << "\n*** Recomputing the Hessian representation\n";
      // fill the D matrix (diagonal)
      this->D_matrix[this->current_memory_slot] = 10.; // TODO dot(s_new, y_new)
      std::cout << "D: "; print_vector(std::cout, this->D_matrix);

      // fill the L matrix (lower triangular with a zero diagonal)
      for (size_t column_index: Range(this->used_memory_size)) {
         for (size_t row_index: Range(column_index+1, this->used_memory_size)) {
            this->L_matrix.entry(row_index, column_index) = static_cast<double>(column_index + row_index); // TODO dot(s_i, y_j)
         }
      }
      std::cout << "L:\n" << this->L_matrix;

      // assemble m x m matrix M = S^T B0 S + L D^{-1} L^T
      // Ltilde = L D^{-1/2}
      DenseMatrix<double> Ltilde_matrix(this->L_matrix); // copy L into Ltilde
      for (size_t column_index: Range(this->used_memory_size)) {
         const double scaling = 1./std::sqrt(this->D_matrix[column_index]);
         for (size_t row_index: Range(column_index+1, this->used_memory_size)) {
            Ltilde_matrix.entry(row_index, column_index) *= scaling;
         }
      }
      std::cout << "Ltilde:\n" << Ltilde_matrix;

      // form M = L D^{-1} L^T + S^T S = Ltilde Ltilde^T + S^T S
      this->M_matrix.clear();
      LBFGSHessian::perform_high_rank_update(this->M_matrix, this->used_memory_size, this->memory_size, Ltilde_matrix,
         this->used_memory_size, this->memory_size);
      LBFGSHessian::perform_high_rank_update_transpose(this->M_matrix, this->used_memory_size, this->memory_size,
         this->S_matrix, this->used_memory_size, this->number_variables);
      std::cout << "M:\n" << this->M_matrix;

      // compute the Cholesky factors J of M = J J^T (J overwrites M)
      LBFGSHessian::compute_cholesky_factors(this->M_matrix, this->used_memory_size, this->memory_size);
      std::cout << "J:\n" << this->M_matrix;
   }

   void LBFGSHessian::perform_high_rank_update(DenseMatrix<double>& matrix, size_t matrix_dimension, size_t matrix_leading_dimension,
         DenseMatrix<double>& high_rank_correction, size_t correction_rank, size_t correction_leading_dimension) {
      std::cout << "Performing rank " << correction_rank << " update\n";
      char uplo = 'L'; // lower triangular
      char trans = 'N';
      int n = static_cast<int>(matrix_dimension); // dimension of matrix
      int k = static_cast<int>(correction_rank); // number of columns of high_rank_correction
      int lda = static_cast<int>(correction_leading_dimension); // number of rows of high_rank_correction
      int ldc = static_cast<int>(matrix_leading_dimension); // leading dimension of matrix
      double alpha = 1.;
      double beta = 1.;
      assert(lda >= std::max(1, n) && "LBFGSHessian::perform_high_rank_update assumption on lda is violated");
      assert(ldc >= std::max(1, n) && "LBFGSHessian::perform_high_rank_update assumption on ldc is violated");
      LAPACK_symmetric_high_rank_update(&uplo, &trans, &n, &k, &alpha, high_rank_correction.data(), &lda, &beta,
         matrix.data(), &ldc);
   }

   void LBFGSHessian::perform_high_rank_update_transpose(DenseMatrix<double>& matrix, size_t matrix_dimension, size_t matrix_leading_dimension,
         DenseMatrix<double>& high_rank_correction, size_t correction_rank, size_t correction_leading_dimension) {
      std::cout << "Performing rank " << correction_rank << " update\n";
      char uplo = 'L'; // lower triangular
      char trans = 'T';
      int n = static_cast<int>(matrix_dimension); // dimension of matrix
      int k = static_cast<int>(correction_leading_dimension); // number of rows of high_rank_correction
      int lda = static_cast<int>(correction_leading_dimension); // number of rows of high_rank_correction
      int ldc = static_cast<int>(matrix_leading_dimension); // leading dimension of matrix
      double alpha = 1.;
      double beta = 1.;
      assert(lda >= std::max(1, k) && "LBFGSHessian::perform_high_rank_update_transpose assumption on lda is violated");
      assert(ldc >= std::max(1, n) && "LBFGSHessian::perform_high_rank_update assumption on ldc is violated");
      LAPACK_symmetric_high_rank_update(&uplo, &trans, &n, &k, &alpha, high_rank_correction.data(), &lda, &beta,
         matrix.data(), &ldc);
   }

   void LBFGSHessian::compute_cholesky_factors(DenseMatrix<double>& matrix, size_t dimension, size_t leading_dimension) {
      char uplo = 'L';
      int info = 0;
      int dimension_int = static_cast<int>(dimension);
      int leading_dimension_int = static_cast<int>(leading_dimension);
      LAPACK_cholesky_factorization(&uplo, &dimension_int, matrix.data(), &leading_dimension_int, &info);
      std::cout << "Cholesky info: " << info << '\n';
      assert(info == 0 && "The Cholesky factorization failed");
   }
} // namespace