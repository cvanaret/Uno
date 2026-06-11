// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include <limits>
#include "LSR1Hessian.hpp"
#include "linear_algebra/LAPACK_extension.hpp"
#include "model/Model.hpp"
#include "options/Options.hpp"
#include "symbolic/Range.hpp"
#include "symbolic/ScalarMultiple.hpp"
#include "symbolic/Triangular.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   LSR1Hessian::LSR1Hessian(const Model& model, double objective_multiplier, const Options& options):
         DirectQuasiNewtonHessian("L-SR1", model, objective_multiplier, options),
         pivot_max_magnitude(options.get_double("LSR1_pivot_max_magnitude")),
         LD(this->memory_size, this->memory_size),
         N(this->memory_size, this->memory_size),
         U(this->model.number_variables, this->memory_size),
         correction_scaling(this->memory_size) {
   }

   bool LSR1Hessian::is_positive_definite() const {
      return false;
   }

   void LSR1Hessian::initialize_statistics(Statistics& statistics) const {
      statistics.add_column("|SR1|", Statistics::double_width - 3, 2);
      statistics.set("|SR1|", this->number_entries_in_memory);
   }

   void LSR1Hessian::notify_trial_iterate(Statistics& statistics, const Iterate& current_iterate, const Iterate& trial_iterate,
         EvaluationCache& evaluation_cache) {
      // compute the candidate pair (s, y) WITHOUT modifying the memory yet
      this->compute_candidate_pair(current_iterate, trial_iterate, evaluation_cache);

      // safeguard: if dot(sk, yk) is too small relative to sk and yk, skip the update
      const double norm_sk = norm_2(this->latest_s);
      const double norm_yk = norm_2(this->latest_y);
      // tolerance is √(machine epsilon)
      if (dot(this->latest_s, this->latest_y) < std::sqrt(std::numeric_limits<double>::epsilon()) * norm_sk * norm_yk) {
         DEBUG << "dot(sk, yk) is too small, skipping the update\n";
      }
      else {
         // the curvature test passed; the candidate is validated (N factorization) and committed lazily in
         // recompute_hessian_representation, before its next use. It may still be rejected there if N is singular.
         this->hessian_recomputation_required = true;
      }
      statistics.set("|SR1|", this->number_entries_in_memory);
   }

   // Hessian-vector product where the Hessian approximation is Bk = B0 + U P⁻¹ Uᵀ and B0 = δ I
   // Bk v = (B0 + U P⁻¹ Uᵀ) v = δ v + U (P⁻¹ (Uᵀ v))
   void LSR1Hessian::compute_hessian_vector_product(View<const double> /*x*/, View<const double> vector,
         double objective_multiplier, const Vector<double>& /*constraint_multipliers*/, View<double> result) {
      if (objective_multiplier != this->fixed_objective_multiplier) {
         throw std::runtime_error("The L-SR1 Hessian model was initialized with a different objective multiplier");
      }

      // a recomputation of the Hessian representation may be required
      if (this->hessian_recomputation_required) {
         this->recompute_hessian_representation();
         this->hessian_recomputation_required = false;
      }

      // diagonal contribution δ I
      for (size_t variable_index: Range(this->model.number_variables)) {
         result[variable_index] = this->delta * vector[variable_index];
      }

      // rank-1 contribution: U in R^{n x m}
      // work on each column of U (P⁻¹ (Uᵀ v))
      DEBUG << "U = " << this->U << '\n';
      for (size_t column_index: Range(this->number_entries_in_memory)) {
         const auto current_U_column = this->U.column(column_index);
         double U_coefficient = dot(current_U_column, vector) / this->get_correction_column_scaling(column_index);
         assert(!std::isnan(U_coefficient));
         // result += coefficient * column(P⁻¹) * current_column
         blas1::add(this->model.number_variables, U_coefficient, current_U_column.data(), result.data());
      }
   }

   size_t LSR1Hessian::get_correction_rank() const {
      return this->number_entries_in_memory;
   }

   View<const double> LSR1Hessian::get_correction_column(size_t column_index) const {
      return this->U.column(column_index);
   }

   double LSR1Hessian::get_correction_column_scaling(size_t column_index) const {
      return this->correction_scaling[column_index];
   }

   // protected member functions

   void LSR1Hessian::shift_memory_entries() {
      // shift the memory entries S and Y (drop the oldest, slot 0)
      QuasiNewtonHessian::shift_memory_entries();
      // shift the symmetric matrix LD = D + L + Lᵀ (stored as lower triangular, diagonal included)
      this->shift_lower_triangle(this->LD, true);
   }

   void LSR1Hessian::recompute_hessian_representation() {
      assert(this->hessian_recomputation_required);
      // A candidate pair (latest_s, latest_y) passed the curvature test in notify_trial_iterate. Tentatively augment the
      // memory with it as the newest entry (dropping the oldest if the memory is full) and check that the inner matrix N
      // factorizes; only then physically commit the candidate. This way a skipped (singular-N) update does not displace
      // the existing memory.
      const bool memory_full = (this->number_entries_in_memory == this->memory_size);
      const size_t provisional_number_entries = memory_full ? this->memory_size : this->number_entries_in_memory + 1;
      const size_t newest = provisional_number_entries - 1;
      // when the memory is full, a kept tentative entry p corresponds to the stored slot p + 1 (slot 0 is dropped)
      const size_t offset = memory_full ? 1 : 0;
      const double sTy = dot(this->latest_s, this->latest_y);

      /* build the symmetric matrix N = (D + L + Lᵀ) - δ Sᵀ S for the candidate-augmented memory, in chronological order,
       * WITHOUT modifying the stored memory. N is workspace: only the lower triangle is filled (read by the factorization),
       * and a failed factorization leaves only N clobbered, not the committed representation (U, delta, correction_scaling).
       * Stored LD is also chronological, so its (i+offset, j+offset) entry is the inner product s_{i+off}ᵀ y_{j+off}. */
      auto Nk = this->N.submatrix(provisional_number_entries, provisional_number_entries);
      // leading block: the existing entries we keep
      for (size_t row_index: Range(newest)) {
         for (size_t column_index: Range(row_index + 1)) { // lower triangle, diagonal included
            const double sy = this->LD.entry(row_index + offset, column_index + offset);
            const double ss = dot(this->S.column(row_index + offset), this->S.column(column_index + offset));
            this->N.entry(row_index, column_index) = sy - this->delta * ss;
         }
      }
      // last row: the candidate against the kept entries, and the candidate diagonal
      for (size_t column_index: Range(newest)) {
         const double sy = dot(this->latest_s, this->Y.column(column_index + offset));
         const double ss = dot(this->latest_s, this->S.column(column_index + offset));
         this->N.entry(newest, column_index) = sy - this->delta * ss;
      }
      this->N.entry(newest, newest) = sTy - this->delta * dot(this->latest_s, this->latest_s);
      DEBUG << "> N: " << this->N;

      /* compute an LDLᵀ (signed Cholesky) factorization without pivoting of N */
      // N = J P Jᵀ with J lower triangular with unit diagonal, and P diagonal, indefinite with nonzero elements
      // J and P overwrite N
      const bool ldlt_success = ldlt_nopiv_lvl2_rightlooking(Nk.data(), Nk.number_rows, Nk.leading_dimension,
         this->pivot_max_magnitude);
      // if the factorization failed with near-0 pivot, skip the update: the stored memory is untouched and U, delta and
      // correction_scaling still describe the last committed representation, so the next product remains correct
      if (!ldlt_success) {
         DEBUG << "N is singular, skipping the update\n";
         return;
      }
      DEBUG << "The matrix N was successfully factorized\n";

      /* the update is valid: physically commit the candidate as the newest entry (shifting out the oldest if full) */
      const size_t committed_slot = this->commit_memory_entry();
      assert(committed_slot == newest);
      // set the new last row of the symmetric LD (lower triangle, diagonal included)
      for (size_t column_index: Range(committed_slot)) {
         this->LD.entry(committed_slot, column_index) = dot(this->S.column(committed_slot), this->Y.column(column_index));
      }
      this->LD.entry(committed_slot, committed_slot) = sTy;
      DEBUG << "> LD: " << this->LD;

      /* update the initial Hessian approximation δ I */
      this->delta = this->compute_delta();
      DEBUG << "Initial identity multiple: " << this->delta << "\n";

      /* form U = (Y - δ S) J⁻ᵀ */
      auto Uk = this->U.submatrix(this->model.number_variables, this->number_entries_in_memory);
      for (size_t column_index: Range(this->number_entries_in_memory)) {
         this->U.column(column_index) = this->Y.column(column_index) - this->delta*this->S.column(column_index);
      }
      Uk *= transpose(inverse(lower_unit_triangular(Nk)));
      DEBUG << "> U: " << this->U;

      /* publish the column scalings P (diagonal of N = J P Jᵀ); only now, on a committed update, is the representation
       * replaced, so a subsequently skipped update leaves these last valid values in place */
      for (size_t column_index: Range(this->number_entries_in_memory)) {
         this->correction_scaling[column_index] = this->N.entry(column_index, column_index);
      }
   }

   double LSR1Hessian::compute_delta() const {
      assert(0 < this->number_entries_in_memory);
      // TODO safeguard
      return 1.;
   }
} // namespace
