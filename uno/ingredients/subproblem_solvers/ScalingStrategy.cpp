// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include "ScalingStrategy.hpp"
#include "COOLinearSystem.hpp"
#include "symbolic/Range.hpp"
#include "tools/Logger.hpp"

namespace uno {
   ScalingStrategy::ScalingStrategy(const LinearSystemScalingParameters& parameters):
         parameters(parameters),
         active(parameters.method != LinearSystemScalingMethod::NONE && !parameters.on_demand) {
   }

   void ScalingStrategy::begin_new_matrix() {
      this->new_matrix = true;
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
      // the scaling factors are computed at most once per matrix
      if (!this->scaling_up_to_date) {
         this->compute_scaling(linear_system);
      }
      const size_t number_nonzeros = linear_system.number_nonzeros;
      double* values = linear_system.matrix_values.data();
      this->unscaled_values.assign(values, values + number_nonzeros);
      for (size_t nonzero_index: Range(number_nonzeros)) {
         values[nonzero_index] *= this->entry_factors[nonzero_index];
      }
      return true;
   }

   void ScalingStrategy::unscale_matrix(COOLinearSystem& linear_system) const {
      std::copy(this->unscaled_values.begin(), this->unscaled_values.end(), linear_system.matrix_values.data());
   }

   void ScalingStrategy::notify_factorization(std::optional<size_t> number_delayed_pivots, size_t dimension) {
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
      this->scaling_method.compute(this->parameters.method, linear_system.matrix_values.data(), this->scaling);
      this->scaling_method.compute_entry_factors(this->scaling, this->entry_factors);
      this->scaling_up_to_date = true;
      ++this->number_scaling_computations;
      if (!this->scaling.empty()) {
         DEBUG << "Linear system scaling computed (range [" << *std::min_element(this->scaling.begin(), this->scaling.end()) <<
            ", " << *std::max_element(this->scaling.begin(), this->scaling.end()) << "])\n";
      }
   }
} // namespace
