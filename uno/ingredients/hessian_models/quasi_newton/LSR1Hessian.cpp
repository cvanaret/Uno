// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "LSR1Hessian.hpp"
#include "linear_algebra/LAPACK_extension.hpp"
#include "model/Model.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   LSR1Hessian::LSR1Hessian(const Model& model, double objective_multiplier, const Options& options):
            QuasiNewtonHessian("L-SR1", model, objective_multiplier, options),
         LD(this->memory_size, this->memory_size),
         N(this->memory_size, this->memory_size),
         U(this->model.number_variables, this->memory_size) {
   }

   bool LSR1Hessian::is_positive_definite() const {
      return false;
   }

   void LSR1Hessian::initialize_statistics(Statistics& statistics) const {
      statistics.add_column("|SR1|", Statistics::double_width - 3, 2, Statistics::column_order.at("|SR1|"));
      statistics.set("|SR1|", this->number_entries_in_memory);
   }

   void LSR1Hessian::notify_accepted_iterate(Statistics& statistics, const Iterate& current_iterate, const Iterate& trial_iterate,
         EvaluationCache& evaluation_cache) {
      statistics.set("|SR1|", this->number_entries_in_memory);
      DEBUG << "\n*** Adding entries to the SR1 memory at slot " << this->current_index << '\n';
      // update the matrices S and Y
      this->update_memory_entries(current_iterate, trial_iterate, evaluation_cache);
      // take the new point into account for the computations
      const size_t size = std::min(this->number_entries_in_memory + 1, this->memory_size);

      /* update the entries of the lower triangular matrix LD := D + L + Lᵀ */
      // the entries of LD depend on this->current_index (1 row and 1 column, possibly empty)
      // row this->current_index
      for (size_t column_index: Range(this->current_index)) {
         this->LD.entry(this->current_index, column_index) = dot(this->S.column(this->current_index), this->Y.column(column_index));
      }
      // column this->current_index (including the diagonal)
      for (size_t row_index: Range(this->current_index, size)) {
         this->LD.entry(row_index, this->current_index) = dot(this->S.column(row_index), this->Y.column(this->current_index));
      }
      DEBUG << "> LD: " << this->LD;

      /* build the symmetric matrix N = D + L + Lᵀ - δ Sᵀ S */
      const auto Sk = this->S.submatrix(this->model.number_variables, size);
      auto Lk = this->LD.submatrix(size, size);
      auto Nk = this->N.submatrix(size, size);
      Nk = (-this->delta) * (transpose(Sk)*Sk);
      for (size_t column_index: Range(size)) {
         for (size_t row_index: Range(column_index, size)) {
            this->N.entry(row_index, column_index) += this->LD.entry(row_index, column_index);
         }
      }
      DEBUG << "> N: " << this->N;

      /* compute a LDL' (signed Cholesky) factorization without pivoting of N */
      // TODO because the factorization may fail, make a copy
      const double near_zero_pivot_tolerance = 0.; // TODO
      const bool ldlt_success = ldlt_nopiv_lvl2_rightlooking(Nk.data(), Nk.number_rows, Nk.leading_dimension, near_zero_pivot_tolerance);
      if (ldlt_success) {
         this->validate_update();
      }
      else {
         DEBUG << "Skipping the update\n";
      }
      statistics.set("|SR1|", this->number_entries_in_memory);
      // notify_accepted_iterate is called at the end of a major iteration. Since we don't know yet whether the
      // Hessian approximation will be used, we delay the update to the beginning of the next major iteration
      this->hessian_recomputation_required = true;
   }

   // Hessian-vector product where the Hessian approximation is Bk = B0 + U P⁻¹ Uᵀ and B0 = δ I
   // Bk v = (B0 + U P⁻¹ Uᵀ) v = δ v + U (P⁻¹ (Uᵀ v))
   void LSR1Hessian::compute_hessian_vector_product(const double* /*x*/, const double* vector,
         double objective_multiplier, const Vector<double>& /*constraint_multipliers*/, double* result) {
      if (objective_multiplier != this->fixed_objective_multiplier) {
         throw std::runtime_error("The L-SR1 Hessian model was initialized with a different objective multiplier");
      }

      // update_limited_memory() has updated the limited memory. A recomputation of the Hessian representation may be required
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
      for (size_t column_index: Range(this->number_entries_in_memory)) {
         const auto current_U_column = this->U.column(column_index);
         double U_coefficient = dot(current_U_column, vector);
         U_coefficient *= 1.; // TODO
         // result += coefficient * column(P⁻¹) * current_column
         blas1::add(this->model.number_variables, U_coefficient, current_U_column.data(), result);
      }
   }

   size_t LSR1Hessian::get_correction_rank() const {
      return this->number_entries_in_memory;
   }

   VectorView<std::vector<double>> LSR1Hessian::get_correction_column(size_t column_index) const {
      return this->U.column(column_index);
   }

   double LSR1Hessian::get_correction_column_scaling(size_t column_index) const {
      return this->N.entry(column_index, column_index);
   }

   // protected member functions

   void LSR1Hessian::recompute_hessian_representation() {
      // include the new entry into account in the computations
      const size_t provisional_number_entries = std::min(this->number_entries_in_memory + 1, this->memory_size);

      /* update the entries of the symmetric matrix LD := D + L + Lᵀ (represented as lower triangular) */
      // the entries of LD depend on this->current_index (1 row and 1 column)
      // row this->current_index
      for (size_t column_index: Range(this->current_index)) {
         this->LD.entry(this->current_index, column_index) = dot(this->S.column(this->current_index), this->Y.column(column_index));
      }
      // column this->current_index (including the diagonal)
      for (size_t row_index: Range(this->current_index, provisional_number_entries)) {
         this->LD.entry(row_index, this->current_index) = dot(this->S.column(row_index), this->Y.column(this->current_index));
      }
      DEBUG << "> LD: " << this->LD;

      /* build the symmetric matrix N = D + L + Lᵀ - δ Sᵀ S */
      const auto Sk = this->S.submatrix(this->model.number_variables, provisional_number_entries);
      auto Nk = this->N.submatrix(provisional_number_entries, provisional_number_entries);
      this->N = this->LD;
      Nk += (-this->delta) * (transpose(Sk)*Sk);
      DEBUG << "> N: " << this->N;

      /* compute a LDL' (signed Cholesky) factorization without pivoting of N */
      // N = J P Jᵀ with J lower triangular with unit diagonal, and P diagonal, indefinite with nonzero elements
      // J and P overwrite N
      constexpr double near_zero_pivot_tolerance = 0.; // TODO
      const bool ldlt_success = ldlt_nopiv_lvl2_rightlooking(Nk.data(), Nk.number_rows, Nk.leading_dimension,
         near_zero_pivot_tolerance);
      // if the factorization failed with near-0 pivot, skip the update
      if (!ldlt_success) {
         DEBUG << "Skipping the update\n";
         return;
      }

      this->validate_update();
      DEBUG << "The matrix N was successfully factorized\n";
      const auto& Jk = Nk; // TODO signal unit diagonal

      /* update the initial Hessian approximation δ I */
      this->delta = this->compute_delta();
      DEBUG << "Initial identity multiple: " << this->delta << "\n";

      /* form U = (Y - δ S) J⁻ᵀ */
      const auto Yk = this->Y.submatrix(this->model.number_variables, this->number_entries_in_memory);
      auto Uk = this->U.submatrix(this->model.number_variables, this->number_entries_in_memory);
      // Uk = Yk - this->delta * Sk;
      // Uk *= transpose(inverse(Jk));

      // increment the slot: if we exceed the size of the memory, we start over and replace the oldest point in memory
      this->current_index = (this->current_index + 1) % this->memory_size;
   }

   // compute δ = yᵀ y / sᵀ y at the last entry
   // (7.20) in Numerical optimization (Nocedal & Wright)
   double LSR1Hessian::compute_delta() const {
      assert(0 < this->number_entries_in_memory);
      const auto current_column_S = this->S.column(this->current_index);
      const auto current_column_Y = this->Y.column(this->current_index);
      const double numerator = dot(current_column_Y, current_column_Y);
      const double denominator = dot(current_column_S, current_column_Y); // TODO get term from this->LD
      return numerator/denominator;
   }
} // namespace