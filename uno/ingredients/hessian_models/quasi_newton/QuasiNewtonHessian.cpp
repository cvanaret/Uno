// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "QuasiNewtonHessian.hpp"
#include "model/Model.hpp"
#include "optimization/EvaluationCache.hpp"
#include "optimization/Iterate.hpp"
#include "options/Options.hpp"
#include "symbolic/Range.hpp"
#include "symbolic/Subtraction.hpp"

namespace uno {
   QuasiNewtonHessian::QuasiNewtonHessian(const std::string_view name, const Model& model, double objective_multiplier,
      const Options& options):
         HessianModel(name),
         model(model),
         fixed_objective_multiplier(objective_multiplier),
         memory_size(options.get_unsigned_int("quasi_newton_memory_size")),
         S(this->model.number_variables, this->memory_size),
         Y(this->model.number_variables, this->memory_size),
         current_lagrangian_gradient(this->model.number_variables),
         trial_lagrangian_gradient(this->model.number_variables),
         latest_s(this->model.number_variables),
         latest_y(this->model.number_variables) {
      if (this->memory_size <= 0) {
         throw std::runtime_error("The quasi-Newton memory size should be positive");
      }
   }

   bool QuasiNewtonHessian::has_hessian_operator() const {
      return true;
   }

   bool QuasiNewtonHessian::has_curvature() const {
      return true;
   }

   // protected member functions

   void QuasiNewtonHessian::compute_candidate_pair(const Iterate& current_iterate, const Iterate& trial_iterate,
         EvaluationCache& evaluation_cache) {
      this->compute_latest_s(current_iterate, trial_iterate);
      this->compute_latest_y(current_iterate, trial_iterate, evaluation_cache);
      DEBUG << "Computed candidate (s, y) pair\n";
      DEBUG << "> s: " << this->latest_s;
      DEBUG << "> y: " << this->latest_y;
   }

   void QuasiNewtonHessian::compute_latest_s(const Iterate& current_iterate, const Iterate& trial_iterate) {
      // TODO check that the S entry isn't too small
      this->latest_s = view(trial_iterate.primals, 0, this->model.number_variables) -
         view(current_iterate.primals, 0, this->model.number_variables);
   }

   // fill the latest y: y = \nabla L(x_k, y_k, z_k) - \nabla L(x_{k-1}, y_k, z_k)
   void QuasiNewtonHessian::compute_latest_y(const Iterate& current_iterate, const Iterate& trial_iterate, EvaluationCache& evaluation_cache) {
      // evaluate Lagrangian gradients at the current and trial iterates, both with the trial multipliers trial_iterate.multipliers
      this->model.evaluate_lagrangian_gradient(current_iterate.primals, trial_iterate.multipliers, this->fixed_objective_multiplier,
         evaluation_cache.current_evaluations, this->current_lagrangian_gradient);
      this->model.evaluate_lagrangian_gradient(trial_iterate.primals, trial_iterate.multipliers, this->fixed_objective_multiplier,
         evaluation_cache.trial_evaluations, this->trial_lagrangian_gradient);
      this->latest_y = this->trial_lagrangian_gradient - this->current_lagrangian_gradient;
   }

   size_t QuasiNewtonHessian::commit_memory_entry() {
      size_t newest_slot;
      if (this->number_entries_in_memory == this->memory_size) {
         // the memory is full: physically drop the oldest entry (slot 0) and append the new one at the end
         this->shift_memory_entries();
         newest_slot = this->memory_size - 1;
      }
      else {
         // the memory has room: append the new entry at the end
         newest_slot = this->number_entries_in_memory;
         ++this->number_entries_in_memory;
      }
      // store the candidate pair in the newest slot
      this->S.column(newest_slot) = this->latest_s;
      this->Y.column(newest_slot) = this->latest_y;
      DEBUG << "Committed (s, y) at slot " << newest_slot << " (" << this->number_entries_in_memory << " entries in memory)\n";
      DEBUG << "> S: " << this->S;
      DEBUG << "> Y: " << this->Y;
      return newest_slot;
   }

   void QuasiNewtonHessian::shift_memory_entries() {
      // only ever called when the memory is full, so all memory_size slots hold valid entries
      // drop the oldest entry (slot 0): slots 1..memory_size-1 move to 0..memory_size-2
      for (size_t slot: Range(this->memory_size - 1)) {
         this->S.column(slot) = this->S.column(slot + 1);
         this->Y.column(slot) = this->Y.column(slot + 1);
      }
   }

   void QuasiNewtonHessian::shift_lower_triangle(DenseMatrix<double>& matrix, bool include_diagonal) const {
      // entry(i, j) = entry(i+1, j+1) over the lower triangle. Reading from higher indices while writing to lower ones
      // (in increasing order) is safe in place.
      for (size_t row_index: Range(this->memory_size - 1)) {
         const size_t column_count = include_diagonal ? row_index + 1 : row_index; // columns [0, column_count)
         for (size_t column_index: Range(column_count)) {
            matrix.entry(row_index, column_index) = matrix.entry(row_index + 1, column_index + 1);
         }
      }
   }

   void QuasiNewtonHessian::shift_upper_triangle(DenseMatrix<double>& matrix, bool include_diagonal) const {
      // entry(i, j) = entry(i+1, j+1) over the upper triangle
      for (size_t row_index: Range(this->memory_size - 1)) {
         const size_t first_column = include_diagonal ? row_index : row_index + 1;
         for (size_t column_index: Range(first_column, this->memory_size - 1)) {
            matrix.entry(row_index, column_index) = matrix.entry(row_index + 1, column_index + 1);
         }
      }
   }
} // namespace
