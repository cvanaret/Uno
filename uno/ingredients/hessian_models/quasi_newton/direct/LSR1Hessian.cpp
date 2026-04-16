// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <limits>
#include "LSR1Hessian.hpp"
#include "linear_algebra/LAPACK_extension.hpp"
#include "model/Model.hpp"
#include "options/Options.hpp"
#include "symbolic/ScalarMultiple.hpp"
#include "symbolic/Triangular.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   LSR1Hessian::LSR1Hessian(const Model& model, double objective_multiplier, const Options& options):
         DirectQuasiNewtonHessian("L-SR1", model, objective_multiplier, options),
         pivot_max_magnitude(options.get_double("LSR1_pivot_max_magnitude")),
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

   void LSR1Hessian::notify_trial_iterate(Statistics& statistics, const Iterate& current_iterate, const Iterate& trial_iterate,
         EvaluationCache& evaluation_cache) {
      statistics.set("|SR1|", this->number_entries_in_memory);
      DEBUG << "\n*** Adding entries to the SR1 memory at slot " << this->current_index << '\n';
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
      DEBUG << "U = " << this->U << '\n';
      for (size_t column_index: Range(this->number_entries_in_memory)) {
         const auto current_U_column = this->U.column(column_index);
         double U_coefficient = dot(current_U_column, vector) / this->get_correction_column_scaling(column_index);
         assert(!std::isnan(U_coefficient));
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

      /* compute an LDL' (signed Cholesky) factorization without pivoting of N */
      // N = J P Jᵀ with J lower triangular with unit diagonal, and P diagonal, indefinite with nonzero elements
      // J and P overwrite N
      const bool ldlt_success = ldlt_nopiv_lvl2_rightlooking(Nk.data(), Nk.number_rows, Nk.leading_dimension,
         this->pivot_max_magnitude);
      // if the factorization failed with near-0 pivot, skip the update
      if (!ldlt_success) {
         DEBUG << "N is singular, skipping the update\n";
         return;
      }

      this->validate_update();
      DEBUG << "The matrix N was successfully factorized\n";

      /* update the initial Hessian approximation δ I */
      this->delta = this->compute_delta();
      DEBUG << "Initial identity multiple: " << this->delta << "\n";

      /* form U = (Y - δ S) J⁻ᵀ */
      const auto Yk = this->Y.submatrix(this->model.number_variables, this->number_entries_in_memory);
      auto Uk = this->U.submatrix(this->model.number_variables, this->number_entries_in_memory);
      for (size_t column_index: Range(this->number_entries_in_memory)) {
         this->U.column(column_index) = this->Y.column(column_index) - this->delta*this->S.column(column_index);
      }
      Uk *= transpose(inverse(lower_unit_triangular(Nk)));

      // increment the slot: if we exceed the size of the memory, we start over and replace the oldest point in memory
      this->current_index = (this->current_index + 1) % this->memory_size;
   }

   double LSR1Hessian::compute_delta() const {
      // we cannot use the L-BFGS initialization because it kills U
      return 1.;
   }
} // namespace