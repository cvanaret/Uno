// Copyright (c) 2025-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "LBFGSHessian.hpp"
#include "model/Model.hpp"
#include "options/Options.hpp"
#include "linear_algebra/BLAS.hpp"
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
   LBFGSHessian::LBFGSHessian(const Model& model, double objective_multiplier, const Options& options):
         DirectQuasiNewtonHessian("L-BFGS", model, objective_multiplier, options),
         // matrices
         L(this->memory_size, this->memory_size),
         D(this->memory_size),
         invsqrt_D(this->memory_size),
         L_invsqrt_D(this->memory_size, this->memory_size),
         M(this->memory_size, this->memory_size),
         U(this->model.number_variables, this->memory_size),
         V(this->model.number_variables, this->memory_size),
         max_skips_before_reset(options.get_unsigned_int("LBFGS_max_skips_before_reset")) {
   }

   bool LBFGSHessian::is_positive_definite() const {
      return true;
   }

   void LBFGSHessian::initialize_statistics(Statistics& statistics) const {
      statistics.add_column("|BFGS|", Statistics::double_width - 2, 2);
      statistics.set("|BFGS|", this->number_entries_in_memory);
   }

   void LBFGSHessian::notify_trial_iterate(Statistics& statistics, const Iterate& current_iterate, const Iterate& trial_iterate,
         Evaluations& current_evaluations, Evaluations& trial_evaluations) {
      // compute the candidate pair (s, y) WITHOUT modifying the memory yet, so that a skipped update is a true no-op
      this->compute_candidate_pair(current_iterate, trial_iterate, current_evaluations, trial_evaluations);

      // safeguard: if dot(sk, yk) is too small relative to sk and yk, skip the update
      // TODO compare against Procedure 18.2 (Damped BFGS Updating) from Numerical Optimization
      const double norm_sk = norm_2(this->latest_s);
      const double norm_yk = norm_2(this->latest_y);
      const double sTy = dot(this->latest_s, this->latest_y);
      // tolerance is √(machine epsilon)
      if (sTy <= std::sqrt(std::numeric_limits<double>::epsilon()) * norm_sk * norm_yk) {
         DEBUG << "dot(sk, yk) is too small, skipping the update\n";
         ++this->consecutive_skips;
         if (this->consecutive_skips >= this->max_skips_before_reset) {
            // reset the limited memory
            DEBUG << "Update was skipped " << this->max_skips_before_reset << " consecutive times, resetting the limited memory\n";
            this->number_entries_in_memory = 0;
            this->consecutive_skips = 0;
            // nothing to recompute: an empty memory yields Bk = δ I
            this->hessian_recomputation_required = false;
         }
         // a skip that does not reset leaves the memory (and therefore the current representation) unchanged
      }
      else {
         DEBUG << "Update is valid\n";
         this->consecutive_skips = 0;
         // append the candidate pair as the newest memory entry. If the memory is full, the oldest entry is physically
         // shifted out (commit_memory_entry -> shift_memory_entries) so that the slots stay in chronological order.
         const size_t newest = this->commit_memory_entry();
         this->D[newest] = sTy;
         // compute the new (last) row of L: L(newest, j) = s_newest · y_j for j < newest. No column update is needed:
         // the newest pair only contributes lower-triangular interactions, which all lie in this last row.
         for (size_t column_index: Range(newest)) {
            this->L.entry(newest, column_index) = dot(this->S.column(newest), this->Y.column(column_index));
         }
         // notify_accepted_iterate is called at the end of a major iteration. Since we don't know yet whether the Hessian
         // approximation will be used, we delay the recomputation of the factored representation to its next use.
         this->hessian_recomputation_required = true;
      }
      statistics.set("|BFGS|", this->number_entries_in_memory);
   }

   // Hessian-vector product where the Hessian approximation is Bk = B0 - U Uᵀ + V Vᵀ and B0 = δ I
   // Bk v = (B0 - U Uᵀ + V Vᵀ) v = δ v - U (Uᵀ v) + V (Vᵀ v)
   void LBFGSHessian::compute_hessian_vector_product(const double* /*x*/, const double* vector,
         double objective_multiplier, const Vector<double>& /*constraint_multipliers*/, double* result) {
      if (objective_multiplier != this->fixed_objective_multiplier) {
         throw std::runtime_error("The quasi-Newton Hessian model was initialized with a different objective multiplier");
      }

      // recompute the Hessian representation if the limited memory was updated
      if (this->hessian_recomputation_required) {
         this->recompute_hessian_representation();
         this->hessian_recomputation_required = false;
      }

      // diagonal contribution δ I
      for (size_t variable_index: Range(this->model.number_variables)) {
         result[variable_index] = this->delta * vector[variable_index];
      }

      // rank-2 contribution
      // (U, V) in R^{n x m}
      // work on each column of U (Uᵀ v)
      for (size_t column_index: Range(this->number_entries_in_memory)) {
         const auto current_U_column = this->U.column(column_index);
         const double U_coefficient = -dot(current_U_column, vector); // minus sign for U
         // result += coefficient * current_column
         blas1::add(this->model.number_variables, U_coefficient, current_U_column.data(), result);
      }
      // work on each column of V (Vᵀ v)
      for (size_t column_index: Range(this->number_entries_in_memory)) {
         const auto current_V_column = this->V.column(column_index);
         const double V_coefficient = dot(current_V_column, vector); // plus sign for V
         // result += coefficient * current_column
         blas1::add(this->model.number_variables, V_coefficient, current_V_column.data(), result);
      }
   }

   size_t LBFGSHessian::get_correction_rank() const {
      return 2 * this->number_entries_in_memory;
   }

   // get a column of (U V)
   VectorView<const double> LBFGSHessian::get_correction_column(size_t column_index) const {
      if (column_index < this->number_entries_in_memory) {
         return this->U.column(column_index);
      }
      else {
         return this->V.column(column_index - this->number_entries_in_memory);
      }
   }

   // get the scaling of a given column (-1 for U, +1 for V)
   double LBFGSHessian::get_correction_column_scaling(size_t column_index) const {
      if (column_index < this->number_entries_in_memory) {
         return -1.;
      }
      else {
         return 1.;
      }
   }

   // protected member functions

   void LBFGSHessian::shift_memory_entries() {
      // shift the memory entries S and Y (drop the oldest, slot 0)
      QuasiNewtonHessian::shift_memory_entries();
      // shift the incrementally maintained cached quantities accordingly: the diagonal D and the strictly lower
      // triangular L. The scalings invsqrt_D, L_invsqrt_D and V are not shifted; they are recomputed from scratch in
      // recompute_hessian_representation.
      for (size_t slot: Range(this->memory_size - 1)) {
         this->D[slot] = this->D[slot + 1];
      }
      this->shift_lower_triangle(this->L, false); // L is strictly lower triangular (no diagonal)
   }

   void LBFGSHessian::recompute_hessian_representation() {
      assert(0 < this->number_entries_in_memory);

      DEBUG << "\n*** Recomputing the Hessian representation with " << this->number_entries_in_memory << " entries\n";
      // some matrices were allocated with the maximum memory size. We will work with submatrices instead (of size
      // this->number_entries_in_memory, not this->memory_size)
      const auto Sk = this->S.submatrix(this->model.number_variables, this->number_entries_in_memory);
      const auto L_invsqrt_Dk = this->L_invsqrt_D.submatrix(this->number_entries_in_memory, this->number_entries_in_memory);
      auto Mk = this->M.submatrix(this->number_entries_in_memory, this->number_entries_in_memory);

      /* update the initial Hessian approximation δ I */
      this->delta = this->compute_delta();
      DEBUG << "Initial identity multiple: " << this->delta << "\n";

      /* recompute the diagonal scaling invsqrt_D = D^{-1/2} */
      for (size_t index: Range(this->number_entries_in_memory)) {
         this->invsqrt_D[index] = 1./std::sqrt(this->D[index]);
      }

      /* recompute L_invsqrt_D = L D^{-1/2} (strictly lower triangular) from the chronological L.
       * L is maintained incrementally (one new row per accepted update, shifted on wraparound); the diagonal and upper
       * triangle stay zero, so M = L_invsqrt_D L_invsqrt_Dᵀ is correct. */
      for (size_t row_index: Range(this->number_entries_in_memory)) {
         for (size_t column_index: Range(row_index)) {
            this->L_invsqrt_D.entry(row_index, column_index) = this->invsqrt_D[column_index] * this->L.entry(row_index, column_index);
         }
      }
      DEBUG << "> L: " << this->L;
      DEBUG << "> L_invsqrt_D: " << this->L_invsqrt_D;

      /* recompute V = Y D^{-1/2} */
      for (size_t index: Range(this->number_entries_in_memory)) {
         this->V.column(index) = this->invsqrt_D[index] * this->Y.column(index);
      }

      /* form M = L D⁻¹ Lᵀ + Sᵀ B0 S = L_invsqrt_D L_invsqrt_Dᵀ + δ Sᵀ S */
      Mk = L_invsqrt_Dk * transpose(L_invsqrt_Dk);
      Mk += this->delta * (transpose(Sk) * Sk);
      DEBUG << "> M: " << this->M;

      /* compute the Cholesky factor J of M = J Jᵀ */
      const bool success = Mk.compute_cholesky_factorization(); // J overwrites M
      DEBUG << "Cholesky success: " << success << '\n';
      DEBUG << "> J: " << this->M;

      /* form U = (δ S + Y D⁻¹ Lᵀ) J⁻ᵀ = (δ S + V L_invsqrt_Dᵀ) J⁻ᵀ */
      const auto Jk = this->M.submatrix(this->number_entries_in_memory, this->number_entries_in_memory); // J overwrites M
      auto Uk = this->U.submatrix(this->model.number_variables, this->number_entries_in_memory);
      const auto Vk = this->V.submatrix(this->model.number_variables, this->number_entries_in_memory);
      Uk = Sk;
      Uk = this->delta * Uk + Vk * transpose(L_invsqrt_Dk);
      Uk *= transpose(inverse(lower_triangular(Jk)));
      DEBUG << "> U: " << this->U;
      DEBUG << "> V: " << this->V << '\n';

      // note: the newest entry already lives at slot (number_entries_in_memory - 1); there is no circular index to advance
   }

   // compute (guarded) δ = sᵀy / sᵀs at the newest (last) entry
   double LBFGSHessian::compute_delta() const {
      assert(0 < this->number_entries_in_memory);
      const size_t newest = this->number_entries_in_memory - 1;
      const double sTy = this->D[newest];
      const auto s = this->S.column(newest);
      const double sTs = dot(s, s);
      return std::max(this->delta_lower_bound, std::min(this->delta_upper_bound, sTy/sTs));
   }
} // namespace
