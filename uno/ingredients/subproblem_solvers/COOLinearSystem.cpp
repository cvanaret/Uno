// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "COOLinearSystem.hpp"
#include "ingredients/subproblem/Subproblem.hpp"

namespace uno {
   COOLinearSystem::COOLinearSystem(int solver_indexing): solver_indexing(solver_indexing) {
   }

   void COOLinearSystem::initialize_hessian(const Subproblem& subproblem) {
      const size_t number_hessian_nonzeros = subproblem.number_hessian_nonzeros();
      const size_t number_primal_inertia_correction_nonzeros = subproblem.number_primal_inertia_correction_nonzeros();
      // Hessian
      this->dimension = subproblem.number_variables;
      this->number_nonzeros = number_hessian_nonzeros + number_primal_inertia_correction_nonzeros;
      this->matrix_row_indices.resize(this->number_nonzeros);
      this->matrix_column_indices.resize(this->number_nonzeros);
      // compute the COO sparse representation
      subproblem.compute_regularized_hessian_sparsity(this->matrix_row_indices.data(), this->matrix_column_indices.data(),
         this->solver_indexing);
      this->matrix_values.resize(this->number_nonzeros);
      this->rhs.resize(this->dimension);
      this->solution.resize(this->dimension);
      // create the views
      this->hessian = view(this->matrix_values.data(), number_hessian_nonzeros);
      this->jacobian = view(this->hessian.end(), 0);
      this->primal_inertia_correction = view(this->jacobian.end(), number_primal_inertia_correction_nonzeros);
      this->dual_inertia_correction = view(this->primal_inertia_correction.end(), 0);
   }

   void COOLinearSystem::initialize_augmented_system(const Subproblem& subproblem) {
      const size_t number_hessian_nonzeros = subproblem.number_hessian_nonzeros();
      const size_t number_jacobian_nonzeros = subproblem.number_jacobian_nonzeros();
      const size_t number_primal_inertia_correction_nonzeros = subproblem.number_primal_inertia_correction_nonzeros();
      const size_t number_dual_inertia_correction_nonzeros = subproblem.number_dual_inertia_correction_nonzeros();
      // augmented system
      this->dimension = subproblem.number_variables + subproblem.number_constraints;
      this->number_nonzeros = number_hessian_nonzeros + number_jacobian_nonzeros + number_primal_inertia_correction_nonzeros +
         number_dual_inertia_correction_nonzeros;
      this->matrix_row_indices.resize(this->number_nonzeros);
      this->matrix_column_indices.resize(this->number_nonzeros);
      // compute the COO sparse representation
      subproblem.compute_regularized_augmented_matrix_sparsity(this->matrix_row_indices.data(), this->matrix_column_indices.data(),
         this->solver_indexing);
      this->matrix_values.resize(this->number_nonzeros);
      this->rhs.resize(this->dimension);
      this->solution.resize(this->dimension);
      // create the views
      this->hessian = view(this->matrix_values.data(), number_hessian_nonzeros);
      this->jacobian = view(this->hessian.end(), number_jacobian_nonzeros);
      this->primal_inertia_correction = view(this->jacobian.end(), number_primal_inertia_correction_nonzeros);
      this->dual_inertia_correction = view(this->primal_inertia_correction.end(), number_dual_inertia_correction_nonzeros);
   }

   double COOLinearSystem::compute_hessian_quadratic_form(const Subproblem& /*subproblem*/,
         const Vector<double>& /*vector*/) const {
      return 0.;
   }
} // namespace
