// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "LBFGSHessian.hpp"
#include "model/Model.hpp"
#include "options/Options.hpp"
#include "linear_algebra/LAPACK.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/Range.hpp"
#include "tools/Statistics.hpp"

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
      std::cout << "Adding vector to L-BFGS memory at slot " << this->current_memory_slot << '\n';
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

      // TODO atm we use delta * I
      for (size_t variable_index: Range(model.number_variables)) {
         result[variable_index] = this->initial_identity_multiple * vector[variable_index];
      }
   }

   std::string LBFGSHessian::get_name() const {
      return "L-BFGS";
   }

   // protected member functions

   void LBFGSHessian::update_memory(const Model& model, Iterate& current_iterate, Iterate& trial_iterate) {
      DEBUG << "*** Updating the BFGS memory\n";
      // update the matrices Y, S and D
      this->update_Y_matrix(model, current_iterate, trial_iterate);
      // TODO check that the S entry isn't too small
      this->update_S_matrix(current_iterate, trial_iterate);
      this->update_D_matrix();
      DEBUG << "> Y: " << this->Y_matrix;
      DEBUG << "> S: " << this->S_matrix;
      DEBUG << "> D: "; print_vector(DEBUG, this->D_matrix);

      // check that the latest D entry s^T y is > 0
      if (0. < this->D_matrix[this->current_memory_slot]) {
         DEBUG << "Adding vector to L-BFGS memory at slot " << this->current_memory_slot << '\n';
         this->number_entries_in_memory = std::min(this->number_entries_in_memory + 1, this->memory_size);
         this->hessian_recomputation_required = true;
         DEBUG << "There are now " << this->number_entries_in_memory << " iterates in memory (capacity " <<
            this->memory_size << ")\n";
      }
      else {
         DEBUG << "Skipping the update\n";
      }
   }

   void LBFGSHessian::update_Y_matrix(const Model& model, Iterate& current_iterate, Iterate& trial_iterate) {
      // fill the Y matrix: y = \nabla L(x_k, y_k, z_k) - \nabla L(x_{k-1}, y_k, z_k)
      model.evaluate_objective_gradient(current_iterate.primals, current_iterate.evaluations.objective_gradient);
      model.evaluate_constraint_jacobian(current_iterate.primals, current_iterate.evaluations.constraint_jacobian);
      model.evaluate_objective_gradient(trial_iterate.primals, trial_iterate.evaluations.objective_gradient);
      model.evaluate_constraint_jacobian(trial_iterate.primals, trial_iterate.evaluations.constraint_jacobian);

      // evaluate Lagrangian gradients at the current and trial iterates, both with the trial multipliers
      LagrangianGradient<double> current_lagrangian_gradient(this->dimension);
      LagrangianGradient<double> trial_lagrangian_gradient(this->dimension);
      model.evaluate_lagrangian_gradient(current_lagrangian_gradient, current_iterate, trial_iterate.multipliers);
      model.evaluate_lagrangian_gradient(trial_lagrangian_gradient, trial_iterate, trial_iterate.multipliers);
      // TODO objective multiplier is hardcoded for the moment
      const auto assembled_current_lagrangian_gradient = current_lagrangian_gradient.assemble(1.);
      const auto assembled_trial_lagrangian_gradient = trial_lagrangian_gradient.assemble(1.);
      for (size_t variable_index: Range(this->dimension)) {
         this->Y_matrix.entry(variable_index, this->current_memory_slot) = assembled_trial_lagrangian_gradient[variable_index] -
            assembled_current_lagrangian_gradient[variable_index];
      }
   }

   void LBFGSHessian::update_S_matrix(const Iterate& current_iterate, const Iterate& trial_iterate) {
      for (size_t variable_index: Range(this->dimension)) {
         this->S_matrix.entry(variable_index, this->current_memory_slot) = trial_iterate.primals[variable_index] -
            current_iterate.primals[variable_index];
      }
   }

   void LBFGSHessian::update_D_matrix() {
      this->D_matrix[this->current_memory_slot] = dot(this->S_matrix.column(this->current_memory_slot),
         this->Y_matrix.column(this->current_memory_slot));
   }

   void LBFGSHessian::recompute_hessian_representation() {
      DEBUG << "\n*** Recomputing the Hessian representation\n";
      // TODO figure out if we're extending or replacing in memory
      // fill the L matrix (lower triangular with a zero diagonal)
      for (size_t column_index: Range(this->number_entries_in_memory)) {
         for (size_t row_index: Range(column_index+1, this->number_entries_in_memory)) {
            this->L_matrix.entry(row_index, column_index) = dot(this->S_matrix.column(row_index), this->Y_matrix.column(column_index));
         }
      }
      DEBUG << "> L: " << this->L_matrix;

      // form Ltilde = L D^{-1/2}
      // TODO preallocate
      DenseMatrix<double> Ltilde_matrix(this->L_matrix); // copy L into Ltilde
      for (size_t column_index: Range(this->number_entries_in_memory)) {
         const double scaling = 1./std::sqrt(this->D_matrix[column_index]);
         for (size_t row_index: Range(column_index+1, this->number_entries_in_memory)) {
            Ltilde_matrix.entry(row_index, column_index) *= scaling;
         }
      }
      DEBUG << "> Ltilde: " << Ltilde_matrix;

      // update the initial Hessian approximation delta * I, where delta = 1/gamma and gamma = s^T y / y^T y at the last entry
      this->initial_identity_multiple = (0 < this->number_entries_in_memory) ? this->compute_initial_identity_factor() : 1.;
      DEBUG << "Initial identity multiple: " << this->initial_identity_multiple << "\n";

      // form m x m matrix M = L D^{-1} L^T + S^T B0 S = Ltilde Ltilde^T + S^T B0 S
      this->M_matrix.clear();
      // add M += Ltilde Ltilde^T
      LBFGSHessian::perform_high_rank_update(this->M_matrix, this->number_entries_in_memory, this->memory_size, Ltilde_matrix,
         this->number_entries_in_memory, this->memory_size, 1., 1.);
      // add M += S^T B0 S (= delta S^T S)
      LBFGSHessian::perform_high_rank_update_transpose(this->M_matrix, this->number_entries_in_memory, this->memory_size,
         this->S_matrix, this->number_entries_in_memory, this->dimension, this->initial_identity_multiple, 1.);
      DEBUG << "> M: " << this->M_matrix;
      // compute the Cholesky factors J of M = J J^T (J overwrites M)
      LBFGSHessian::compute_cholesky_factors(this->M_matrix, this->number_entries_in_memory, this->memory_size);
      DEBUG << "> J: " << this->M_matrix << '\n';

      // increment the slot: if we exceed the size of the memory, we start over and replace the older point in memory
      this->current_memory_slot = (this->current_memory_slot + 1) % this->memory_size;
   }

   // precondition: 0 < this->number_entries_in_memory
   double LBFGSHessian::compute_initial_identity_factor() const {
      // TODO get the correct columns
      const auto last_column_Y = this->Y_matrix.column(this->current_memory_slot);
      const auto last_column_S = this->S_matrix.column(this->current_memory_slot);
      // return delta = 1/gamma where gamma is given by (7.20) in Numerical optimization (Nocedal & Wright)
      const double numerator = dot(last_column_Y, last_column_Y);
      const double denominator = dot(last_column_S, last_column_Y); // TODO should be the current D entry
      assert(denominator != 0 && "LBFGSHessian::compute_initial_identity_factor: the denominator is 0");
      return numerator/denominator;
   }

   // performs symmetric rank k update
   // C = alpha A A^T + beta C
   void LBFGSHessian::perform_high_rank_update(DenseMatrix<double>& matrix, size_t matrix_dimension, size_t matrix_leading_dimension,
         DenseMatrix<double>& high_rank_correction, size_t correction_rank, size_t correction_leading_dimension, double alpha, double beta) {
      DEBUG << "Performing rank " << correction_rank << " update\n";
      char uplo = 'L'; // lower triangular
      char trans = 'N';
      int n = static_cast<int>(matrix_dimension); // dimension of matrix
      int k = static_cast<int>(correction_rank); // number of columns of high_rank_correction
      int lda = static_cast<int>(correction_leading_dimension); // number of rows of high_rank_correction
      int ldc = static_cast<int>(matrix_leading_dimension); // leading dimension of matrix
      assert(lda >= std::max(1, n) && "LBFGSHessian::perform_high_rank_update assumption on lda is violated");
      assert(ldc >= std::max(1, n) && "LBFGSHessian::perform_high_rank_update assumption on ldc is violated");
      LAPACK_symmetric_high_rank_update(&uplo, &trans, &n, &k, &alpha, high_rank_correction.data(), &lda, &beta,
         matrix.data(), &ldc);
   }

   // performs symmetric rank k update
   // C = alpha A^T A + beta C
   void LBFGSHessian::perform_high_rank_update_transpose(DenseMatrix<double>& matrix, size_t matrix_dimension, size_t matrix_leading_dimension,
         DenseMatrix<double>& high_rank_correction, size_t correction_rank, size_t correction_leading_dimension, double alpha, double beta) {
      DEBUG << "Performing rank " << correction_rank << " update\n";
      char uplo = 'L'; // lower triangular
      char trans = 'T';
      int n = static_cast<int>(matrix_dimension); // dimension of matrix
      int k = static_cast<int>(correction_leading_dimension); // number of rows of high_rank_correction
      int lda = static_cast<int>(correction_leading_dimension); // number of rows of high_rank_correction
      int ldc = static_cast<int>(matrix_leading_dimension); // leading dimension of matrix
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
      DEBUG << "Cholesky info: " << info << '\n';
   }
} // namespace