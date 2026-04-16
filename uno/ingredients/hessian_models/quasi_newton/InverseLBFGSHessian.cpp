// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "InverseLBFGSHessian.hpp"
#include "model/Model.hpp"
#include "options/Options.hpp"
#include "linear_algebra/BLAS.hpp"
#include "linear_algebra/Vector.hpp"
#include "linear_algebra/VectorView.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/Inverse.hpp"
#include "symbolic/Multiplication.hpp"
#include "symbolic/Range.hpp"
#include "symbolic/ScalarMultiple.hpp"
#include "symbolic/Sum.hpp"
#include "symbolic/Transpose.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   InverseLBFGSHessian::InverseLBFGSHessian(const Model& model, double objective_multiplier, const Options& options):
         QuasiNewtonHessian("inverse L-BFGS", model, objective_multiplier, options),
         // matrices
         L(this->memory_size, this->memory_size),
         D(this->memory_size),
         invsqrt_D(this->memory_size),
         L_invsqrt_D(this->memory_size, this->memory_size),
         M(this->memory_size, this->memory_size),
         U(this->model.number_variables, this->memory_size),
         V(this->model.number_variables, this->memory_size),
         delta_upper_bound(options.get_double("LBFGS_delta_upper_bound")) {
      if (this->memory_size <= 0) {
         throw std::runtime_error("The quasi-Newton memory size should be positive");
      }
   }

   bool InverseLBFGSHessian::has_hessian_matrix() const {
      return false;
   }

   size_t InverseLBFGSHessian::number_nonzeros() const {
      throw std::runtime_error("This member function should not be called.");
   }

   void InverseLBFGSHessian::compute_sparsity(uno_int* /*row_indices*/, uno_int* /*column_indices*/, uno_int /*solver_indexing*/) const {
      throw std::runtime_error("This member function should not be called.");
   }

   bool InverseLBFGSHessian::is_positive_definite() const {
      return true;
   }

   void InverseLBFGSHessian::initialize_statistics(Statistics& statistics) const {
      statistics.add_column("|BFGS|", Statistics::double_width - 2, 2, Statistics::column_order.at("|BFGS|"));
      statistics.set("|BFGS|", this->number_entries_in_memory);
   }

   void InverseLBFGSHessian::notify_trial_iterate(Statistics& statistics, const Iterate& current_iterate, const Iterate& trial_iterate,
         EvaluationCache& evaluation_cache) {
      statistics.set("|BFGS|", this->number_entries_in_memory);
      DEBUG << "\n*** Adding entries to the BFGS memory at slot " << this->current_index << '\n';
      // update the matrices S and Y
      this->update_memory_entries(current_iterate, trial_iterate, evaluation_cache);

      // safeguard: if dot(sk, yk) is too small relative to sk and yk, skip the update
      const auto sk = this->S.column(this->current_index);
      const auto yk = this->Y.column(this->current_index);
      const double norm_sk = dot(sk, sk);
      const double norm_yk = dot(yk, yk);
      // tolerance is √(machine epsilon)
      if (dot(sk, yk) < std::sqrt(std::numeric_limits<double>::epsilon()) * norm_sk * norm_yk) {
         DEBUG << "dot(sk, yk) is too small, skipping the update\n";
      }
      else {
         // notify_accepted_iterate is called at the end of a major iteration. Since we don't know yet whether the
         // Hessian approximation will be used, we delay the update to the beginning of the next major iteration
         this->hessian_recomputation_required = true;
      }
   }

   void InverseLBFGSHessian::evaluate_hessian(Statistics& /*statistics*/, const Vector<double>& /*primal_variables*/,
         double /*objective_multiplier*/, const Vector<double>& /*constraint_multipliers*/, double* /*hessian_values*/) {
      throw std::runtime_error("This member function should not be called.");
   }

   void InverseLBFGSHessian::compute_hessian_vector_product(const double* /*x*/, const double* /*vector*/,
         double /*objective_multiplier*/, const Vector<double>& /*constraint_multipliers*/, double* /*result*/) {
      throw std::runtime_error("This member function should not be called.");
   }

   void InverseLBFGSHessian::compute_inverse_hessian_vector_product(const double* /*x*/, const double* vector, double* result) {
      // recompute the Hessian representation if the limited memory was updated
      if (this->hessian_recomputation_required) {
         this->recompute_hessian_representation();
         hessian_recomputation_required = false;
      }

      std::cout << "INVERSE HESSIAN VECTOR PROD with delta = " << this->delta << '\n';
      for (size_t variable_index: Range(this->model.number_variables)) {
         result[variable_index] += this->delta * vector[variable_index];
      }
   }

   // protected member functions

   void InverseLBFGSHessian::update_D() {
      this->D[this->current_index] = dot(this->S.column(this->current_index), this->Y.column(this->current_index));
      DEBUG << "> diag(D): "; print_vector(DEBUG, this->D);
   }

   void InverseLBFGSHessian::recompute_hessian_representation() {
      // check that the latest D entry sᵀ y is > 0
      // TODO implement Procedure 18.2 (Damped BFGS Updating) from Numerical Optimization
      this->update_D();
      if (this->D[this->current_index] <= 0.) {
         DEBUG << "Skipping the update\n";
         return;
      }

      this->validate_update();
      assert(0 < this->number_entries_in_memory);

      DEBUG << "\n*** Recomputing the Hessian representation with " << this->number_entries_in_memory << " entries\n";
      // note: some matrices were allocated with the maximum memory size. We will work with submatrices instead (of size
      // this->number_entries_in_memory, not this->memory_size)
      const auto Sk = this->S.submatrix(this->model.number_variables, this->number_entries_in_memory);
      const auto L_invsqrt_Dk = this->L_invsqrt_D.submatrix(this->number_entries_in_memory, this->number_entries_in_memory);
      auto Mk = this->M.submatrix(this->number_entries_in_memory, this->number_entries_in_memory);

      /* update the initial Hessian approximation δ I */
      this->delta = this->compute_delta();
      DEBUG << "Initial identity multiple: " << this->delta << "\n";

      /* update the entries of L and L_invsqrt_D = L D^{-1/2} */
      // update invsqrt_D = 1/sqrt_D
      this->invsqrt_D[this->current_index] = 1./std::sqrt(this->D[this->current_index]);
      // the entries of L and L_invsqrt_D depend on this->current_index (1 row and 1 column, possibly empty)
      // row this->current_index
      for (size_t column_index: Range(this->current_index)) {
         this->L.entry(this->current_index, column_index) = dot(this->S.column(this->current_index), this->Y.column(column_index));
         this->L_invsqrt_D.entry(this->current_index, column_index) = this->invsqrt_D[column_index] * this->L.entry(this->current_index, column_index);
      }
      // column this->current_index (excluding the diagonal)
      for (size_t row_index: Range(this->current_index+1, this->number_entries_in_memory)) {
         this->L.entry(row_index, this->current_index) = dot(this->S.column(row_index), this->Y.column(this->current_index));
         this->L_invsqrt_D.entry(row_index, this->current_index) = this->invsqrt_D[this->current_index] * this->L.entry(row_index, this->current_index);
      }
      DEBUG << "> L: " << this->L;
      DEBUG << "> L_invsqrt_D: " << this->L_invsqrt_D;

      /* form M = L D⁻¹ Lᵀ + Sᵀ B0 S = L_invsqrt_D L_invsqrt_Dᵀ + δ Sᵀ S */
      Mk = L_invsqrt_Dk * transpose(L_invsqrt_Dk);
      Mk += this->delta * (transpose(Sk) * Sk);
      DEBUG << "> M: " << this->M;

      /* compute the Cholesky factor J of M = J Jᵀ */
      const bool success = Mk.compute_cholesky_factorization(); // J overwrites M
      DEBUG << "Cholesky success: " << success << '\n';
      DEBUG << "> J: " << this->M;

      /* form V and U */
      // update the current column of V = Y D^{-1/2}
      this->V.column(this->current_index) = this->invsqrt_D[this->current_index] * this->Y.column(this->current_index);
      // form U = (δ S + Y D⁻¹ Lᵀ) J⁻ᵀ
      const auto Jk = this->M.submatrix(this->number_entries_in_memory, this->number_entries_in_memory); // J overwrites M
      auto Uk = this->U.submatrix(this->model.number_variables, this->number_entries_in_memory);
      const auto Vk = this->V.submatrix(this->model.number_variables, this->number_entries_in_memory);
      Uk = Sk;
      Uk = this->delta * Uk + Vk * transpose(L_invsqrt_Dk);
      Uk *= transpose(inverse(Jk));
      DEBUG << "> U: " << this->U;
      DEBUG << "> V: " << this->V << '\n';

      // increment the slot: if we exceed the size of the memory, we start over and replace the oldest point in memory
      this->current_index = (this->current_index + 1) % this->memory_size;
   }

   // compute δ = yᵀ y / sᵀ y at the last entry
   // (7.20) in Numerical optimization (Nocedal & Wright)
   double InverseLBFGSHessian::compute_delta() const {
      assert(0 < this->number_entries_in_memory);
      const auto last_column_Y = this->Y.column(this->current_index);
      const double numerator = dot(last_column_Y, last_column_Y);
      const double denominator = this->D[this->current_index]; // > 0 by the update rule
      return denominator/numerator;
   }
} // namespace