// Copyright (c) 2025-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "LBFGSHessian.hpp"
#include "model/Model.hpp"
#include "options/Options.hpp"
#include "linear_algebra/BLAS.hpp"
#include "linear_algebra/LAPACK.hpp"
#include "linear_algebra/Vector.hpp"
#include "linear_algebra/VectorView.hpp"
#include "optimization/EvaluationCache.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/Range.hpp"
#include "symbolic/ScalarMultiple.hpp"
#include "symbolic/Subtraction.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   LBFGSHessian::LBFGSHessian(const Model& model, double objective_multiplier, const Options& options):
         HessianModel("L-BFGS"),
         model(model),
         fixed_objective_multiplier(objective_multiplier),
         memory_size(options.get_unsigned_int("quasi_newton_memory_size")),
         // matrices
         S(this->model.number_variables, this->memory_size),
         Y(this->model.number_variables, this->memory_size),
         L(this->memory_size, this->memory_size),
         Ltilde(this->memory_size, this->memory_size),
         D(this->memory_size),
         invsqrt_D(this->memory_size),
         M(this->memory_size, this->memory_size),
         U(this->model.number_variables, this->memory_size),
         V(this->model.number_variables, this->memory_size),
         current_lagrangian_gradient(this->model.number_variables),
         trial_lagrangian_gradient(this->model.number_variables) {
      if (this->memory_size <= 0) {
         throw std::runtime_error("The quasi-Newton memory size should be positive");
      }
   }

   bool LBFGSHessian::has_hessian_operator() const {
      return true;
   }

   bool LBFGSHessian::has_hessian_matrix() const {
      // never form the explicit, dense matrix
      return false;
   }

   bool LBFGSHessian::has_curvature() const {
      return true;
   }

   size_t LBFGSHessian::number_nonzeros() const {
      throw std::runtime_error("LBFGSHessian::number_nonzeros should not be called");
   }

   void LBFGSHessian::compute_sparsity(int* /*row_indices*/, int* /*column_indices*/, int /*solver_indexing*/) const {
      throw std::runtime_error("LBFGSHessian::compute_sparsity should not be called");
   }

   bool LBFGSHessian::is_positive_definite() const {
      return true;
   }

   void LBFGSHessian::initialize_statistics(Statistics& statistics) const {
      statistics.add_column("|BFGS|", Statistics::double_width - 2, 2, Statistics::column_order.at("|BFGS|"));
      statistics.set("|BFGS|", this->number_entries_in_memory);
   }

   void LBFGSHessian::notify_accepted_iterate(Statistics& statistics, const Iterate& current_iterate, const Iterate& trial_iterate,
         EvaluationCache& evaluation_cache) {
      DEBUG << "\n*** Updating the BFGS memory at slot " << this->current_index << '\n';
      // update the matrices S, Y and D
      this->update_S(current_iterate, trial_iterate);
      this->update_Y(current_iterate, trial_iterate, evaluation_cache);
      this->update_D();
      DEBUG << "> S: " << this->S;
      DEBUG << "> Y: " << this->Y;
      DEBUG << "> diag(D): "; print_vector(DEBUG, this->D);

      // check that the latest D entry s^T y is > 0
      // TODO implement Procedure 18.2 (Damped BFGS Updating) from Numerical Optimization
      if (0. < this->D[this->current_index]) {
         DEBUG << "S, Y and D updated at slot " << this->current_index << '\n';
         this->number_entries_in_memory = std::min(this->number_entries_in_memory + 1, this->memory_size);
         // notify_accepted_iterate is called at the end of a major iteration. At this point, we don't know yet whether
         // the L-BFGS Hessian will be used. We therefore delay the update to the beginning of the next major iteration
         // by setting the flag hessian_recomputation_required
         this->hessian_recomputation_required = true;
         DEBUG << "There are now " << this->number_entries_in_memory << " entries in memory (capacity " << this->memory_size << ")\n";
      }
      else {
         DEBUG << "Skipping the update\n";
      }
      statistics.set("|BFGS|", this->number_entries_in_memory);
   }

   void LBFGSHessian::evaluate_hessian(Statistics& /*statistics*/, const Vector<double>& /*primal_variables*/,
         double /*objective_multiplier*/, const Vector<double>& /*constraint_multipliers*/, double* /*hessian_values*/) {
      throw std::runtime_error("LBFGSHessian::evaluate_hessian should not be called");
   }

   // Hessian-vector product where the Hessian approximation is Bk = B0 - U U^T + V V^T and B0 = delta I
   // Bk v = (B0 - U U^T + V V^T) v = delta v - U (U^T v) + V (V^T v)
   void LBFGSHessian::compute_hessian_vector_product(const double* /*x*/, const double* vector,
         double objective_multiplier, const Vector<double>& /*constraint_multipliers*/, double* result) {
      if (objective_multiplier != this->fixed_objective_multiplier) {
         throw std::runtime_error("The L-BFGS Hessian model was initialized with a different objective multiplier");
      }

      // update_limited_memory() has updated the limited memory. A recomputation of the Hessian representation may be required
      if (this->hessian_recomputation_required) {
         this->recompute_hessian_representation();
         this->hessian_recomputation_required = false;
      }

      // diagonal contribution delta*I
      for (size_t variable_index: Range(model.number_variables)) {
         result[variable_index] = this->delta * vector[variable_index];
      }

      // rank-2 contribution
      // (U, V) in R^{n x m}
      int n = static_cast<int>(this->model.number_variables);
      int increment = 1;
      // U U^T v
      for (size_t column_index: Range(this->number_entries_in_memory)) {
         const auto current_U_column = this->U.column(column_index);
         const double U_coefficient = -dot(current_U_column, vector); // minus sign for U
         // result += coefficient * current_column
         BLAS_add_vector(&n, &U_coefficient, current_U_column.data(), &increment, result, &increment);
      }
      // V V^T v
      for (size_t column_index: Range(this->number_entries_in_memory)) {
         const auto current_V_column = this->V.column(column_index);
         const double V_coefficient = dot(current_V_column, vector); // plus sign for V
         // result += coefficient * current_column
         BLAS_add_vector(&n, &V_coefficient, current_V_column.data(), &increment, result, &increment);
      }
   }

   // protected member functions

   void LBFGSHessian::update_S(const Iterate& current_iterate, const Iterate& trial_iterate) {
      // TODO check that the S entry isn't too small
      this->S.column(this->current_index) = view(trial_iterate.primals, 0, this->model.number_variables) -
         view(current_iterate.primals, 0, this->model.number_variables);
   }
   
   // fill the Y matrix: y = \nabla L(x_k, y_k, z_k) - \nabla L(x_{k-1}, y_k, z_k)
   void LBFGSHessian::update_Y(const Iterate& current_iterate, const Iterate& trial_iterate, EvaluationCache& evaluation_cache) {
      // evaluate Lagrangian gradients at the current and trial iterates, both with the trial multipliers trial_iterate.multipliers
      this->model.evaluate_lagrangian_gradient(current_iterate.primals, trial_iterate.multipliers, this->fixed_objective_multiplier,
         evaluation_cache.current_evaluations, this->current_lagrangian_gradient);
      this->model.evaluate_lagrangian_gradient(trial_iterate.primals, trial_iterate.multipliers, this->fixed_objective_multiplier,
         evaluation_cache.trial_evaluations, this->trial_lagrangian_gradient);
      this->Y.column(this->current_index) = this->trial_lagrangian_gradient - this->current_lagrangian_gradient;
   }

   void LBFGSHessian::update_D() {
      this->D[this->current_index] = dot(this->S.column(this->current_index), this->Y.column(this->current_index));
   }

   void LBFGSHessian::recompute_hessian_representation() {
      assert(0 < this->number_entries_in_memory);
      assert(0 < this->D[this->current_index]);

      DEBUG << "\n*** Recomputing the Hessian representation with " << this->number_entries_in_memory << " entries\n";

      // update the initial Hessian approximation delta * I, where delta = 1/gamma and gamma = s^T y / y^T y at the last entry
      this->delta = this->compute_delta();
      DEBUG << "Initial identity multiple: " << this->delta << "\n";

      // update invsqrt_D = 1/sqrt_D
      this->invsqrt_D[this->current_index] = 1./std::sqrt(this->D[this->current_index]);

      // update the entries of L and Ltilde = L D^{-1/2} that depend on this->current_index (1 row and 1 column, possibly empty)
      // row this->current_index
      for (size_t column_index: Range(this->current_index)) {
         this->L.entry(this->current_index, column_index) = dot(this->S.column(this->current_index), this->Y.column(column_index));
         this->Ltilde.entry(this->current_index, column_index) = this->invsqrt_D[column_index] * this->L.entry(this->current_index, column_index);
      }
      // column this->current_index
      for (size_t row_index: Range(this->current_index+1, this->number_entries_in_memory)) {
         this->L.entry(row_index, this->current_index) = dot(this->S.column(row_index), this->Y.column(this->current_index));
         this->Ltilde.entry(row_index, this->current_index) = this->invsqrt_D[this->current_index] * this->L.entry(row_index, this->current_index);
      }
      DEBUG << "> L: " << this->L;
      DEBUG << "> Ltilde: " << this->Ltilde;

      // form m x m matrix M = L D^{-1} L^T + S^T B0 S = Ltilde Ltilde^T + delta S^T S
      LBFGSHessian::perform_high_rank_update(this->M, this->number_entries_in_memory, this->memory_size, this->Ltilde,
         this->number_entries_in_memory, this->memory_size, 1., 0.); // Ltilde Ltilde^T
      LBFGSHessian::perform_high_rank_update_transpose(this->M, this->number_entries_in_memory, this->memory_size,
         this->S, this->number_entries_in_memory, this->model.number_variables, this->delta, 1.); // delta S^T S
      DEBUG << "> M: " << this->M;

      // compute the Cholesky factor J of M = J J^T (overwrites M with J)
      LBFGSHessian::compute_cholesky_factors(this->M, this->number_entries_in_memory, this->memory_size);
      DenseMatrix<double>& J_matrix = this->M;
      DEBUG << "> J: " << J_matrix;

      // update the current column of V = Y D^{-1/2}
      for (size_t row_index: Range(this->model.number_variables)) {
         this->V.entry(row_index, this->current_index) = this->invsqrt_D[this->current_index] * this->Y.entry(row_index, this->current_index);
      }

      // compute V * Ltilde^T in W matrix (B A^T with A = Ltilde, B = V, W stored in U)
      this->U = this->V;
      {
         char side = 'R'; //  B := alpha B op(A)
         char uplo = 'L'; // Ltilde is lower triangular
         char transa = 'T'; // op(A) = A^T
         char diag = 'N';
         int m = static_cast<int>(this->model.number_variables); // number of rows of V
         int n = static_cast<int>(this->number_entries_in_memory); // number of columns of V
         double alpha = 1.;
         int lda = static_cast<int>(this->memory_size); // leading dimension of Ltilde
         int ldb = static_cast<int>(this->model.number_variables); // leading dimension of V
         BLAS_triangular_matrix_matrix_product(&side, &uplo, &transa, &diag, &m, &n, &alpha, this->Ltilde.data(), &lda,
            this->U.data(), &ldb);
      }
      // add delta S to U
      for (size_t column_index: Range(this->number_entries_in_memory)) {
         this->U.column(column_index) += this->delta * this->S.column(column_index);
      }
      DEBUG << "> W: " << this->U;

      // solve U J^T = W wrt U (X op(A) = alpha B with A = J, op(A) = A^T, alpha = 1, B = W)
      // U overwrites W
      {
         char side = 'R'; // X op(A) = alpha B
         char uplo = 'L'; // J is lower triangular
         char transa = 'T'; // op(A) = A^T
         char diag = 'N';
         int m = static_cast<int>(this->model.number_variables); // number of rows of W
         int n = static_cast<int>(this->number_entries_in_memory); // number of columns of W
         double alpha = 1.;
         int lda = static_cast<int>(this->memory_size); // leading dimension of J
         int ldb = static_cast<int>(this->model.number_variables); // leading dimension of W
         BLAS_triangular_back_solve(&side, &uplo, &transa, &diag, &m, &n, &alpha, J_matrix.data(), &lda,
            this->U.data(), &ldb);
      }
      DEBUG << "> U: " << this->U;
      DEBUG << "> V: " << this->V << '\n';

      // increment the slot: if we exceed the size of the memory, we start over and replace the oldest point in memory
      this->current_index = (this->current_index + 1) % this->memory_size;
   }

   // compute delta = 1/gamma where gamma is given by (7.20) in Numerical optimization (Nocedal & Wright)
   // precondition: 0 < this->number_entries_in_memory
   double LBFGSHessian::compute_delta() const {
      assert(0 < this->number_entries_in_memory);
      const auto last_column_Y = this->Y.column(this->current_index);
      const double numerator = dot(last_column_Y, last_column_Y);
      const double denominator = this->D[this->current_index]; // > 0 by the update rule
      return numerator/denominator;
   }

   // performs symmetric rank k update
   // C = alpha A A^T + beta C
   void LBFGSHessian::perform_high_rank_update(DenseMatrix<double>& matrix, size_t matrix_dimension, size_t matrix_leading_dimension,
         DenseMatrix<double, MatrixShape::LOWER_TRIANGULAR>& high_rank_correction, size_t correction_rank, size_t correction_leading_dimension, double alpha, double beta) {
      DEBUG << "Performing rank " << correction_rank << " update\n";
      char uplo = 'L'; // lower triangular
      char trans = 'N';
      int n = static_cast<int>(matrix_dimension); // dimension of matrix
      int k = static_cast<int>(correction_rank); // number of columns of high_rank_correction
      int lda = static_cast<int>(correction_leading_dimension); // number of rows of high_rank_correction
      int ldc = static_cast<int>(matrix_leading_dimension); // leading dimension of matrix
      assert(lda >= std::max(1, n));
      assert(ldc >= std::max(1, n));
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
      assert(lda >= std::max(1, k));
      assert(ldc >= std::max(1, n));
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