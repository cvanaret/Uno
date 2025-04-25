// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "LBFGSHessian.hpp"
#include "model/Model.hpp"
#include "options/Options.hpp"
#include "linear_algebra/LAPACK.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "model/Model.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/Range.hpp"

namespace uno {
   LBFGSHessian::LBFGSHessian(const Options& options): HessianModel(),
         memory_size(options.get_unsigned_int("quasi_newton_memory_size")) {
   }

   bool LBFGSHessian::has_implicit_representation() const {
      return true;
   }

   bool LBFGSHessian::has_explicit_representation() const {
      return false;
   }

   bool LBFGSHessian::has_curvature(const Model& /*model*/) const {
      return true;
   }

   size_t LBFGSHessian::number_nonzeros(const Model& /*model*/) const {
      throw std::runtime_error("LBFGSHessian::number_nonzeros should not be called");
   }

   bool LBFGSHessian::is_positive_definite() const {
      return true;
   }

   void LBFGSHessian::initialize(const Model& model) {
      this->dimension = model.number_variables;
      this->S_matrix = DenseMatrix<double>(this->dimension, this->memory_size);
      this->Y_matrix = DenseMatrix<double>(this->dimension, this->memory_size);
      this->L_matrix = DenseMatrix<double>(this->memory_size, this->memory_size);
      this->D_matrix.resize(this->memory_size);
      this->M_matrix = DenseMatrix<double>(this->memory_size, this->memory_size);
   }

   void LBFGSHessian::initialize_statistics(Statistics& statistics, const Options& options) const {
   }

   void LBFGSHessian::notify_accepted_iterate(const Model& model, Iterate& current_iterate, Iterate& trial_iterate) {
      std::cout << "Adding vector to L-BFGS memory at slot " << this->current_available_slot << '\n';
      // this->current_available_slot lives in [0, this->memory_size)
      this->update_memory(model, current_iterate, trial_iterate);
   }
   
   void LBFGSHessian::evaluate_hessian(Statistics& /*statistics*/, const Model& /*model*/, const Vector<double>& /*primal_variables*/,
         double /*objective_multiplier*/, const Vector<double>& /*constraint_multipliers*/, SymmetricMatrix<size_t, double>& /*hessian*/) {
      throw std::runtime_error("LBFGSHessian::evaluate_hessian should not be called");
   }

   void LBFGSHessian::compute_hessian_vector_product(const Model& model, const double* vector, double /*objective_multiplier*/,
         const Vector<double>& /*constraint_multipliers*/, double* result) {
      if (this->hessian_recomputation_required) {
         this->recompute_hessian_representation();
         this->hessian_recomputation_required = false;
      }

      // TODO atm we use the identity
      for (size_t variable_index: Range(model.number_variables)) {
         result[variable_index] = vector[variable_index];
      }
   }

   std::string LBFGSHessian::get_name() const {
      return "L-BFGS";
   }

   // protected member functions

   void LBFGSHessian::update_memory(const Model& model, Iterate& current_iterate, Iterate& trial_iterate) {
      std::cout << "\n*** Adding vector to L-BFGS memory at slot " << this->current_available_slot << '\n';
      // TODO figure out if we're extending or replacing in memory

      // update the matrices
      // TODO check that the S entry isn't too small
      LBFGSHessian::update_S_matrix(current_iterate, trial_iterate);
      LBFGSHessian::update_Y_matrix(model, current_iterate, trial_iterate);

      // update the bookkeeping
      this->number_iterates_in_memory = std::min(this->number_iterates_in_memory + 1, this->memory_size);
      this->hessian_recomputation_required = true;
      std::cout << "There are now " << this->number_iterates_in_memory << " iterates in memory (capacity " <<
         this->memory_size << ")\n";
   }

   void LBFGSHessian::update_S_matrix(const Iterate& current_iterate, const Iterate& trial_iterate) {
      for (size_t variable_index: Range(this->dimension)) {
         this->S_matrix.entry(variable_index, this->current_available_slot) = trial_iterate.primals[variable_index] -
            current_iterate.primals[variable_index];
      }
      std::cout << "S:\n" << this->S_matrix;
   }

   void LBFGSHessian::update_Y_matrix(const Model& model, Iterate& current_iterate, Iterate& trial_iterate) {
      // fill the Y matrix: y = \nabla L(x_k, y_k, z_k) - \nabla L(x_{k-1}, y_k, z_k)
      model.evaluate_objective_gradient(current_iterate.primals, current_iterate.evaluations.objective_gradient);
      model.evaluate_constraint_jacobian(current_iterate.primals, current_iterate.evaluations.constraint_jacobian);
      model.evaluate_objective_gradient(trial_iterate.primals, trial_iterate.evaluations.objective_gradient);
      model.evaluate_constraint_jacobian(trial_iterate.primals, trial_iterate.evaluations.constraint_jacobian);
   }

   void LBFGSHessian::recompute_hessian_representation() {
      std::cout << "\n*** Recomputing the Hessian representation\n";
      // fill the D matrix (diagonal)
      this->D_matrix[this->current_available_slot] = 10.; // TODO dot(s_new, y_new)
      std::cout << "D: "; print_vector(std::cout, this->D_matrix);

      // fill the L matrix (lower triangular with a zero diagonal)
      for (size_t column_index: Range(this->number_iterates_in_memory)) {
         for (size_t row_index: Range(column_index+1, this->number_iterates_in_memory)) {
            this->L_matrix.entry(row_index, column_index) = 1.; // TODO dot(s_i, y_j)
         }
      }
      std::cout << "L:\n" << this->L_matrix;

      // assemble m x m matrix M = S^T B0 S + L D^{-1} L^T
      // Ltilde = L D^{-1/2}
      DenseMatrix<double> Ltilde_matrix(this->L_matrix); // copy L into Ltilde
      for (size_t column_index: Range(this->number_iterates_in_memory)) {
         const double scaling = 1./std::sqrt(this->D_matrix[column_index]);
         for (size_t row_index: Range(column_index+1, this->number_iterates_in_memory)) {
            Ltilde_matrix.entry(row_index, column_index) *= scaling;
         }
      }
      std::cout << "Ltilde:\n" << Ltilde_matrix;

      // form M = L D^{-1} L^T + S^T S = Ltilde Ltilde^T + S^T S
      this->M_matrix.clear();
      // add M += Ltilde Ltilde^T
      LBFGSHessian::perform_high_rank_update(this->M_matrix, this->number_iterates_in_memory, this->memory_size, Ltilde_matrix,
         this->number_iterates_in_memory, this->memory_size);
      // add M += S^T S
      LBFGSHessian::perform_high_rank_update_transpose(this->M_matrix, this->number_iterates_in_memory, this->memory_size,
         this->S_matrix, this->number_iterates_in_memory, this->dimension);
      std::cout << "M:\n" << this->M_matrix;

      // compute the Cholesky factors J of M = J J^T (J overwrites M)
      LBFGSHessian::compute_cholesky_factors(this->M_matrix, this->number_iterates_in_memory, this->memory_size);
      std::cout << "J:\n" << this->M_matrix;

      // increment the slot: if we exceed the size of the memory, we start over and replace the older point in memory
      this->current_available_slot = (this->current_available_slot + 1) % this->memory_size;
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
      int M_dimension = static_cast<int>(dimension);
      int M_leading_dimension = static_cast<int>(leading_dimension);
      LAPACK_cholesky_factorization(&uplo, &M_dimension, matrix.data(), &M_leading_dimension, &info);
      std::cout << "Cholesky info: " << info << '\n';
   }
} // namespace