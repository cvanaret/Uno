// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "COOLinearSystem.hpp"
#include "ingredients/subproblem/Subproblem.hpp"

namespace uno {
   COOLinearSystem::COOLinearSystem(int solver_indexing): solver_indexing(solver_indexing) {
   }

   void COOLinearSystem::initialize_hessian(const Subproblem& subproblem) {
      // Hessian
      this->dimension = subproblem.number_variables;
      this->number_nonzeros = subproblem.number_regularized_hessian_nonzeros();
      this->matrix_row_indices.resize(this->number_nonzeros);
      this->matrix_column_indices.resize(this->number_nonzeros);
      // compute the COO sparse representation
      subproblem.compute_regularized_hessian_sparsity(this->matrix_row_indices.data(), this->matrix_column_indices.data(),
         this->solver_indexing);
      this->matrix_values.resize(this->number_nonzeros);
      this->rhs.resize(this->dimension);
      this->solution.resize(this->dimension);
   }

   void COOLinearSystem::initialize_augmented_system(const Subproblem& subproblem) {
      // augmented system
      this->dimension = subproblem.number_variables + subproblem.number_constraints;
      this->number_nonzeros = subproblem.number_regularized_augmented_system_nonzeros();
      this->matrix_row_indices.resize(this->number_nonzeros);
      this->matrix_column_indices.resize(this->number_nonzeros);
      // compute the COO sparse representation
      subproblem.compute_regularized_augmented_matrix_sparsity(this->matrix_row_indices.data(), this->matrix_column_indices.data(),
         this->solver_indexing);
      this->matrix_values.resize(this->number_nonzeros);
      this->rhs.resize(this->dimension);
      this->solution.resize(this->dimension);
   }

   double COOLinearSystem::compute_hessian_quadratic_form(const Subproblem& /*subproblem*/,
         const Vector<double>& /*vector*/) const {
      return 0.;
   }
} // namespace
