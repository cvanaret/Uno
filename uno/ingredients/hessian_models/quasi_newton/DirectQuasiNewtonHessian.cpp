// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "DirectQuasiNewtonHessian.hpp"
#include "model/Model.hpp"
#include "options/Options.hpp"

namespace uno {
   DirectQuasiNewtonHessian::DirectQuasiNewtonHessian(const std::string_view name, const Model& model, double objective_multiplier,
      const Options& options):
         QuasiNewtonHessian(name, model, objective_multiplier, options) {
   }

   bool DirectQuasiNewtonHessian::has_hessian_operator() const {
      return true;
   }

   // this will NOT be called by WoodburyEQPSolver
   bool DirectQuasiNewtonHessian::has_hessian_matrix() const {
      // never form the explicit, dense matrix
      return false;
   }

   bool DirectQuasiNewtonHessian::has_curvature() const {
      return true;
   }

   // this can only be called by WoodburyEQPSolver
   size_t DirectQuasiNewtonHessian::number_nonzeros() const {
      // count only the diagonal contribution
      return this->model.number_variables;
   }

   // this can only be called by WoodburyEQPSolver
   void DirectQuasiNewtonHessian::compute_sparsity(int* row_indices, int* column_indices, int solver_indexing) const {
      // diagonal contribution
      for (size_t variable_index: Range(this->model.number_variables)) {
         row_indices[variable_index] = static_cast<int>(variable_index) + solver_indexing;
         column_indices[variable_index] = static_cast<int>(variable_index) + solver_indexing;
      }
   }

   // forms the diagonal part of the L-BFGS Hessian approximation
   // this can only be called by WoodburyEQPSolver
   void DirectQuasiNewtonHessian::evaluate_hessian(Statistics& /*statistics*/, const Vector<double>& /*primal_variables*/,
         double /*objective_multiplier*/, const Vector<double>& /*constraint_multipliers*/, double* hessian_values) {
      // recompute the Hessian representation if the limited memory was updated
      if (this->hessian_recomputation_required) {
         this->recompute_hessian_representation();
         this->hessian_recomputation_required = false;
      }

      // diagonal contribution
      DEBUG << "Setting diagonal contribution of L-BFGS Hessian\n";
      for (size_t variable_index: Range(this->model.number_variables)) {
         hessian_values[variable_index] = this->delta;
      }
   }
} // namespace