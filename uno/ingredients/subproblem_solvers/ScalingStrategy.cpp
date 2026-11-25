// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <algorithm>
#include "ScalingStrategy.hpp"
#include "COOLinearSystem.hpp"
#include "symbolic/Range.hpp"
#include "tools/Logger.hpp"

namespace uno {
   ScalingStrategy::ScalingStrategy(const LinearSystemScalingParameters& parameters):
         parameters(parameters),
         method(parameters.method),
         active(parameters.method != LinearSystemScalingMethod::NONE && !parameters.on_demand) {
   }

   void ScalingStrategy::begin_new_matrix(bool monitor_fill_growth) {
      this->new_matrix = true;
      this->monitor_fill_of_new_matrix = monitor_fill_growth;
   }

   bool ScalingStrategy::activate() {
      if (this->active || this->parameters.method == LinearSystemScalingMethod::NONE) {
         return false;
      }
      if (!this->pending_activation) {
         INFO << "The linear systems will be scaled\n";
      }
      this->active = true;
      this->pending_activation = false;
      this->scaling_up_to_date = false; // computed at the next factorization
      return true;
   }

   bool ScalingStrategy::scale_matrix(COOLinearSystem& linear_system) {
      this->first_factorization_of_matrix = this->new_matrix;
      this->new_matrix = false;
      if (this->first_factorization_of_matrix) {
         this->fill_of_matrix_monitored = this->monitor_fill_of_new_matrix;
      }
      if (this->first_factorization_of_matrix) {
         this->scaling_up_to_date = false;
         if (this->pending_activation) {
            this->pending_activation = false;
            this->active = true;
         }
      }
      this->last_factorization_scaled = this->active;
      if (!this->active) {
         return false;
      }
      // the scaling factors are computed once per matrix, and reused by its refactorizations as long as they keep the
      // scaled entries bounded (see is_scaling_valid_for)
      const size_t number_nonzeros = linear_system.number_nonzeros;
      double* values = linear_system.matrix_values.data();
      if (!this->scaling_up_to_date || !this->is_scaling_valid_for(values, number_nonzeros)) {
         this->compute_scaling(linear_system);
      }
      this->unscaled_values.assign(values, values + number_nonzeros);
      for (size_t nonzero_index: Range(number_nonzeros)) {
         values[nonzero_index] *= this->entry_factors[nonzero_index];
      }
      return true;
   }

   // the scalings bound the entries of the scaled matrix they were computed for (|S K S| <= 1 for the matching-based
   // scaling, ~1 for Ruiz). A refactorization of the same matrix with a regularization (inertia correction) changes its
   // diagonal: since the factors may span tens of orders of magnitude (the barrier terms at the end of interior-point
   // methods), the scaled regularization delta s_i^2 may become huge on some indices and the scaled matrix is then
   // badly scaled (on A0NNDNSL, the fill of the retries grew from 3.4e5 to more than 2.4e7 entries). The scaling is
   // then recomputed. One pass over the entries, negligible compared with a factorization
   bool ScalingStrategy::is_scaling_valid_for(const double* values, size_t number_nonzeros) const {
      for (size_t nonzero_index: Range(number_nonzeros)) {
         if (ScalingStrategy::maximum_scaled_entry_for_reuse < std::abs(values[nonzero_index] * this->entry_factors[nonzero_index])) {
            return false;
         }
      }
      return true;
   }

   void ScalingStrategy::unscale_matrix(COOLinearSystem& linear_system) const {
      std::copy(this->unscaled_values.begin(), this->unscaled_values.end(), linear_system.matrix_values.data());
   }

   void ScalingStrategy::notify_factorization(std::optional<size_t> number_delayed_pivots, size_t dimension,
         std::optional<size_t> number_factor_entries) {
      // fill growth: only the first factorization of each monitored matrix (the refactorizations of the same matrix, e.g.
      // with a regularization, are not representative)
      if (number_factor_entries.has_value() && this->first_factorization_of_matrix && this->fill_of_matrix_monitored) {
         this->monitor_fill_growth(*number_factor_entries);
      }
      // too many delayed pivots in the first (unscaled) factorization of the matrix: switch the scaling on from the next
      // matrix on
      if (this->active || !this->first_factorization_of_matrix || this->parameters.method == LinearSystemScalingMethod::NONE ||
            1. <= this->parameters.delayed_pivots_threshold || !number_delayed_pivots.has_value()) {
         return;
      }
      if (this->parameters.delayed_pivots_threshold * static_cast<double>(dimension) < static_cast<double>(*number_delayed_pivots)) {
         INFO << "The factorization delayed " << *number_delayed_pivots << " pivots (dimension " << dimension <<
            "): the linear systems will be scaled\n";
         this->pending_activation = true;
      }
   }

   // the factors (with at least fill_growth_minimum_entries entries) became much denser than in the sparsest first
   // factorization of a matrix so far: switch on the matching-based scaling from the next factorization on (the solves
   // with the current factorization keep the current scaling). A threshold <= 1 disables the trigger
   void ScalingStrategy::monitor_fill_growth(size_t number_factor_entries) {
      if (this->parameters.method == LinearSystemScalingMethod::NONE || this->parameters.fill_growth_threshold <= 1. ||
            (this->active && this->method == LinearSystemScalingMethod::MATCHING) || number_factor_entries == 0) {
         return; // disabled, already switched, or no meaningful reference (empty factors)
      }
      if (!this->sparsest_factor_entries.has_value() || number_factor_entries < *this->sparsest_factor_entries) {
         this->sparsest_factor_entries = number_factor_entries;
      }
      // the reference is kept for the factors below the minimum number of entries: the fill growth of a matrix that
      // becomes large enough is measured against the sparsest matrices
      else if (this->parameters.fill_growth_minimum_entries <= number_factor_entries &&
            this->parameters.fill_growth_threshold * static_cast<double>(*this->sparsest_factor_entries) <
            static_cast<double>(number_factor_entries)) {
         INFO << "The factors have " << number_factor_entries << " entries (" << *this->sparsest_factor_entries <<
            " at best): the linear systems will be scaled with a matching\n";
         this->method = LinearSystemScalingMethod::MATCHING;
         this->active = true;
         this->pending_activation = false;
         this->scaling_up_to_date = false; // computed at the next factorization
      }
   }

   const double* ScalingStrategy::scale_rhs(const double* rhs, size_t dimension, size_t number_of_rhs) {
      this->scaled_rhs.resize(dimension * number_of_rhs);
      for (size_t column_index: Range(number_of_rhs)) {
         for (size_t index: Range(dimension)) {
            this->scaled_rhs[column_index * dimension + index] = this->scaling[index] * rhs[column_index * dimension + index];
         }
      }
      return this->scaled_rhs.data();
   }

   void ScalingStrategy::scale_solution(double* solution, size_t dimension, size_t number_of_rhs) const {
      for (size_t column_index: Range(number_of_rhs)) {
         for (size_t index: Range(dimension)) {
            solution[column_index * dimension + index] *= this->scaling[index];
         }
      }
   }

   void ScalingStrategy::compute_scaling(const COOLinearSystem& linear_system) {
      // the sparsity pattern is compressed lazily, when the scaling is first needed (again if the pattern changed)
      if (!this->scaling_method.is_initialized() || !this->scaling_method.matches(linear_system.dimension, linear_system.number_nonzeros)) {
         this->scaling_method.initialize(linear_system.dimension, linear_system.matrix_row_indices.data(),
            linear_system.matrix_column_indices.data(), linear_system.number_nonzeros, linear_system.get_solver_indexing());
      }
      this->scaling_method.compute(this->method, linear_system.matrix_values.data(), this->scaling);
      this->scaling_method.compute_entry_factors(this->scaling, this->entry_factors);
      this->scaling_up_to_date = true;
      ++this->number_scaling_computations;
      if (!this->scaling.empty()) {
         DEBUG << "Linear system scaling computed (range [" << *std::min_element(this->scaling.begin(), this->scaling.end()) <<
            ", " << *std::max_element(this->scaling.begin(), this->scaling.end()) << "])\n";
      }
   }
} // namespace
