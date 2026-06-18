// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "InverseLBFGSHessian.hpp"
#include "model/Model.hpp"
#include "options/Options.hpp"
#include "linear_algebra/Vector.hpp"
#include "linear_algebra/VectorView.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/Inverse.hpp"
#include "symbolic/Multiplication.hpp"
#include "symbolic/Range.hpp"
#include "symbolic/ScalarMultiple.hpp"
#include "symbolic/Transpose.hpp"
#include "symbolic/UnaryNegation.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   InverseLBFGSHessian::InverseLBFGSHessian(const Model& model, const Options& options):
         QuasiNewtonHessian("inverse L-BFGS", model, 1., options),
         R(this->memory_size, this->memory_size),
         STv(this->memory_size),
         YTv(this->memory_size) {
   }

   bool InverseLBFGSHessian::has_hessian_matrix() const {
      return false;
   }

   size_t InverseLBFGSHessian::number_nonzeros() const {
      throw std::runtime_error("InverseLBFGSHessian::number_nonzeros should not be called.");
   }

   void InverseLBFGSHessian::compute_sparsity(uno_int* /*row_indices*/, uno_int* /*column_indices*/, uno_int /*solver_indexing*/) const {
      throw std::runtime_error("InverseLBFGSHessian::compute_sparsity should not be called.");
   }

   bool InverseLBFGSHessian::is_positive_definite() const {
      return true;
   }

   void InverseLBFGSHessian::initialize_statistics(Statistics& statistics) const {
      statistics.add_column("|BFGS|", Statistics::double_width - 2, 2);
      statistics.set("|BFGS|", this->number_entries_in_memory);
   }

   void InverseLBFGSHessian::notify_trial_iterate(Statistics& statistics, const Iterate& current_iterate, const Iterate& trial_iterate,
         EvaluationCache& evaluation_cache) {
      // compute the candidate pair (s, y) WITHOUT modifying the memory yet
      this->compute_candidate_pair(current_iterate, trial_iterate, evaluation_cache);

      // safeguard: if dot(sk, yk) is too small relative to sk and yk, skip the update
      const double norm_sk = norm_2(this->latest_s);
      const double norm_yk = norm_2(this->latest_y);
      const double sTy = dot(this->latest_s, this->latest_y);
      // tolerance is √(machine epsilon)
      if (sTy <= std::sqrt(std::numeric_limits<double>::epsilon()) * norm_sk * norm_yk) {
         DEBUG << "dot(sk, yk) is too small, skipping the update\n";
      }
      else {
         // append the candidate pair as the newest memory entry, physically shifting out the oldest if the memory is full
         const size_t newest = this->commit_memory_entry();
         // R is upper triangular with R(i, j) = sᵢ · yⱼ for i ≤ j. The newest pair (last slot) fills the last column,
         // including the diagonal: R(i, newest) = sᵢ · y_newest for i ≤ newest.
         for (size_t row_index: Range(newest)) {
            this->R.entry(row_index, newest) = dot(this->S.column(row_index), this->Y.column(newest));
         }
         this->R.entry(newest, newest) = sTy;
         // notify_accepted_iterate is called at the end of a major iteration. Since we don't know yet whether the
         // Hessian approximation will be used, we delay the recomputation of δ to its next use.
         this->hessian_recomputation_required = true;
      }
      statistics.set("|BFGS|", this->number_entries_in_memory);
   }

   void InverseLBFGSHessian::evaluate_hessian(Statistics& /*statistics*/, const Vector<double>& /*primal_variables*/,
         double /*objective_multiplier*/, const Vector<double>& /*constraint_multipliers*/, double* /*hessian_values*/) {
      throw std::runtime_error("InverseLBFGSHessian::evaluate_hessian should not be called.");
   }

   void InverseLBFGSHessian::compute_hessian_vector_product(const double* /*x*/, const double* /*vector*/,
         double /*objective_multiplier*/, const Vector<double>& /*constraint_multipliers*/, double* /*result*/) {
      throw std::runtime_error("InverseLBFGSHessian::compute_hessian_vector_product should not be called.");
   }

   void InverseLBFGSHessian::compute_inverse_hessian_vector_product(const double* /*x*/, const double* vector, double* result) {
      // recompute the Hessian representation if the limited memory was updated
      if (this->hessian_recomputation_required) {
         this->recompute_hessian_representation();
         hessian_recomputation_required = false;
      }

      const auto v = view(vector, 0, this->model.number_variables);
      auto Hv = view(result, 0, this->model.number_variables);

      // contribution of low-rank corrections
      if (0 < this->number_entries_in_memory) {
         const auto Sk = this->S.submatrix(this->model.number_variables, this->number_entries_in_memory);
         const auto Yk = this->Y.submatrix(this->model.number_variables, this->number_entries_in_memory);
         auto SkTv = view(this->STv, 0, this->number_entries_in_memory);
         auto YkTv = view(this->YTv, 0, this->number_entries_in_memory);
         const auto Rk = this->R.submatrix(this->number_entries_in_memory, this->number_entries_in_memory);

         // compute Sᵀv and Yᵀv
         SkTv = transpose(Sk) * v;
         YkTv = transpose(Yk) * v;
         DEBUG2 << "Sᵀv = "; print_vector(DEBUG2, SkTv);
         DEBUG2 << "Yᵀv = "; print_vector(DEBUG2, YkTv);

         // compute R⁻ᵀ(Yᵀv)
         Vector<double> a(this->number_entries_in_memory);
         a = transpose(inverse(upper_triangular(Rk))) * YkTv;
         DEBUG2 << "R⁻ᵀ(Yᵀv) = " << a << "\n";
         // compute R⁻¹(Sᵀv)
         Vector<double> b(this->number_entries_in_memory);
         b = inverse(upper_triangular(Rk)) * SkTv;
         DEBUG2 << "R⁻¹(Sᵀv) = " << b << "\n";
         Vector<double> p1(this->number_entries_in_memory);
         Vector<double> p2(this->number_entries_in_memory);
         p1 = -this->delta * a;
         p2 = -b;
         Vector<double> c(this->model.number_variables);
         c = Yk * b;
         Vector<double> d(this->number_entries_in_memory);
         d = transpose(Yk) * c;
         d.scale(this->delta);
         // scale b with diag(R)
         for (size_t index: Range(this->number_entries_in_memory)) {
            b[index] *= this->R.entry(index, index);
         }
         d += b;
         Vector<double> e(this->number_entries_in_memory);
         e = transpose(inverse(upper_triangular(Rk))) * d;
         p1 += e;
         // Hv = δv + S p1 + δ Y p2
         Hv = Yk * p2;
         Hv.scale(this->delta);
         Hv += Sk * p1;
      }
      else {
         Hv.fill(0.);
      }

      // diagonal contribution δ I
      Hv += this->delta * v;
   }

   // protected member functions

   void InverseLBFGSHessian::shift_memory_entries() {
      // shift the memory entries S and Y (drop the oldest, slot 0)
      QuasiNewtonHessian::shift_memory_entries();
      // shift the upper triangular R = upper(Sᵀ Y) (diagonal included)
      this->shift_upper_triangle(this->R, true);
   }

   void InverseLBFGSHessian::recompute_hessian_representation() {
      assert(0 < this->number_entries_in_memory);
      DEBUG << "\n*** Recomputing the Hessian representation with " << this->number_entries_in_memory << " entries\n";

      // R is fully maintained in notify_trial_iterate (new last column on each accepted update, shifted on wraparound),
      // so only the initial inverse-Hessian scaling δ needs to be refreshed here
      this->delta = this->compute_delta();
      DEBUG << "Initial identity multiple: " << this->delta << "\n";
      DEBUG << "> R: " << this->R;
   }

   // compute (guarded) δ = sᵀy / yᵀy at the newest (last) entry
   double InverseLBFGSHessian::compute_delta() const {
      assert(0 < this->number_entries_in_memory);
      const size_t newest = this->number_entries_in_memory - 1;
      const double sTy = this->R.entry(newest, newest);
      const auto y = this->Y.column(newest);
      const double yTy = dot(y, y);
      return std::max(this->delta_lower_bound, std::min(this->delta_upper_bound, sTy/yTy));
   }
} // namespace
