// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "QuasiNewtonHessian.hpp"
#include "model/Model.hpp"
#include "optimization/EvaluationCache.hpp"
#include "optimization/Iterate.hpp"
#include "options/Options.hpp"
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
         trial_lagrangian_gradient(this->model.number_variables) {
      if (this->memory_size <= 0) {
         throw std::runtime_error("The quasi-Newton memory size should be positive");
      }
   }

   bool QuasiNewtonHessian::has_hessian_operator() const {
      return true;
   }

   // this will NOT be called by WoodburyEQPSolver
   bool QuasiNewtonHessian::has_hessian_matrix() const {
      // never form the explicit, dense matrix
      return false;
   }

   bool QuasiNewtonHessian::has_curvature() const {
      return true;
   }

   // this can only be called by WoodburyEQPSolver
   size_t QuasiNewtonHessian::number_nonzeros() const {
      // count only the diagonal contribution
      return this->model.number_variables;
   }

   // this can only be called by WoodburyEQPSolver
   void QuasiNewtonHessian::compute_sparsity(int* row_indices, int* column_indices, int solver_indexing) const {
      // diagonal contribution
      for (size_t variable_index: Range(this->model.number_variables)) {
         row_indices[variable_index] = static_cast<int>(variable_index) + solver_indexing;
         column_indices[variable_index] = static_cast<int>(variable_index) + solver_indexing;
      }
   }

   // protected member functions

   void QuasiNewtonHessian::update_memory_entries(const Iterate& current_iterate, const Iterate& trial_iterate,
         EvaluationCache& evaluation_cache) {
      this->update_S(current_iterate, trial_iterate);
      this->update_Y(current_iterate, trial_iterate, evaluation_cache);
      DEBUG << "> S: " << this->S;
      DEBUG << "> Y: " << this->Y;
   }

   void QuasiNewtonHessian::update_S(const Iterate& current_iterate, const Iterate& trial_iterate) {
      // TODO check that the S entry isn't too small
      this->S.column(this->current_index) = view(trial_iterate.primals, 0, this->model.number_variables) -
         view(current_iterate.primals, 0, this->model.number_variables);
   }

   // fill the Y matrix: y = \nabla L(x_k, y_k, z_k) - \nabla L(x_{k-1}, y_k, z_k)
   void QuasiNewtonHessian::update_Y(const Iterate& current_iterate, const Iterate& trial_iterate, EvaluationCache& evaluation_cache) {
      // evaluate Lagrangian gradients at the current and trial iterates, both with the trial multipliers trial_iterate.multipliers
      this->model.evaluate_lagrangian_gradient(current_iterate.primals, trial_iterate.multipliers, this->fixed_objective_multiplier,
         evaluation_cache.current_evaluations, this->current_lagrangian_gradient);
      this->model.evaluate_lagrangian_gradient(trial_iterate.primals, trial_iterate.multipliers, this->fixed_objective_multiplier,
         evaluation_cache.trial_evaluations, this->trial_lagrangian_gradient);
      this->Y.column(this->current_index) = this->trial_lagrangian_gradient - this->current_lagrangian_gradient;
   }

   void QuasiNewtonHessian::validate_update() {
      DEBUG << "S and Y successfully updated at slot " << this->current_index << '\n';
      this->number_entries_in_memory = std::min(this->number_entries_in_memory + 1, this->memory_size);
      // notify_accepted_iterate is called at the end of a major iteration. Since we don't know yet whether the
      // Hessian approximation will be used, we delay the update to the beginning of the next major iteration
      this->hessian_recomputation_required = true;
      DEBUG << "There are now " << this->number_entries_in_memory << " entries in memory (capacity " << this->memory_size << ")\n";
   }
} // namespace