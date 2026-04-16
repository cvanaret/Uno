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

      const auto v = view(vector, 0, this->model.number_variables);

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
         DEBUG2 << "STv = "; print_vector(DEBUG2, SkTv);
         DEBUG2 << "YTv = "; print_vector(DEBUG2, YkTv);

         // compute R⁻ᵀ(Yᵀv)
         Vector<double> a(this->number_entries_in_memory);
         a = transpose(inverse(upper_triangular(Rk))) * YkTv;
         DEBUG2 << "R⁻ᵀ(Yᵀv) = " << a << "\n";
         // compute R⁻¹(Sᵀv)
         Vector<double> b(this->number_entries_in_memory);
         b = inverse(upper_triangular(Rk)) * SkTv;
         DEBUG2 << "R⁻¹(Sᵀv) = " << b << "\n";
         // compute p
         Vector<double> p(2 * this->number_entries_in_memory);
         view(p, 0, this->number_entries_in_memory) = -this->delta * a;
         view(p, this->number_entries_in_memory, 2*this->number_entries_in_memory) = -b;
         //
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
         view(p, 0, this->number_entries_in_memory) += e;
         DEBUG2 << "p = " << p << '\n';
         // result = Hv = δv + S p(1, m) + δ Y p(m, 2m)
         view(result, 0, this->model.number_variables) = Yk * view(p, this->number_entries_in_memory, 2*this->number_entries_in_memory);
         view(result, 0, this->model.number_variables).scale(this->delta);
         view(result, 0, this->model.number_variables) += Sk * view(p, 0, this->number_entries_in_memory);
      }

      // diagonal contribution δ I
      view(result, 0, this->model.number_variables) += this->delta * v;
   }

   // protected member functions

   void InverseLBFGSHessian::recompute_hessian_representation() {
      // check that sᵀy is > 0
      // TODO implement Procedure 18.2 (Damped BFGS Updating) from Numerical Optimization
      const double sTy = dot(this->S.column(this->current_index), this->Y.column(this->current_index));
      if (sTy <= 0.) {
         DEBUG << "Skipping the update\n";
         return;
      }
      this->R.entry(this->current_index, this->current_index) = sTy;

      this->validate_update();
      assert(0 < this->number_entries_in_memory);

      DEBUG << "\n*** Recomputing the Hessian representation with " << this->number_entries_in_memory << " entries\n";

      /* update the initial Hessian approximation δ I */
      this->delta = this->compute_delta();
      DEBUG << "Initial identity multiple: " << this->delta << "\n";

      /* update the entries of R: they depend on this->current_index (1 row and 1 column) */
      // column this->current_index (excluding the diagonal)
      for (size_t row_index: Range(this->current_index)) {
         this->R.entry(row_index, this->current_index) = dot(this->S.column(row_index), this->Y.column(this->current_index));
      }
      // row this->current_index
      for (size_t column_index: Range(this->current_index+1, this->number_entries_in_memory)) {
         this->R.entry(this->current_index, column_index) = dot(this->S.column(this->current_index), this->Y.column(column_index));
      }
      DEBUG << "> R: " << this->R;

      // increment the slot: if we exceed the size of the memory, we start over and replace the oldest point in memory
      this->current_index = (this->current_index + 1) % this->memory_size;
   }

   // compute δ = sᵀ y / yᵀ y at the last entry
   double InverseLBFGSHessian::compute_delta() const {
      assert(0 < this->number_entries_in_memory);
      const auto y = this->Y.column(this->current_index);
      const double sTy = this->R.entry(this->current_index, this->current_index);
      const double yTy = dot(y, y);
      // TODO safeguard
      return sTy/yTy;
   }
} // namespace