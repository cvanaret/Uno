// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include "COOLinearSystem.hpp"
#include "ingredients/subproblem/NewtonSystem.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "symbolic/Range.hpp"

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
      subproblem.compute_regularized_hessian_sparsity(view(this->matrix_row_indices), view(this->matrix_column_indices),
         this->solver_indexing);
      this->matrix_values.resize(this->number_nonzeros);
      this->rhs.resize(this->dimension);
      this->solution.resize(this->dimension);
   }

   void COOLinearSystem::initialize_augmented_system(const Subproblem& subproblem) {
      const size_t number_hessian_nonzeros = subproblem.number_hessian_nonzeros();
      const size_t number_primal_inertia_correction_nonzeros = subproblem.number_variables; // full block
      const size_t number_jacobian_nonzeros = subproblem.number_jacobian_nonzeros();
      const size_t number_dual_inertia_correction_nonzeros = subproblem.number_dual_inertia_correction_nonzeros();
      // augmented system
      this->dimension = subproblem.number_variables + subproblem.number_constraints;
      this->number_nonzeros = number_hessian_nonzeros + number_primal_inertia_correction_nonzeros + number_jacobian_nonzeros
         + number_dual_inertia_correction_nonzeros;
      this->matrix_row_indices.resize(this->number_nonzeros);
      this->matrix_column_indices.resize(this->number_nonzeros);
      // compute the COO sparse representation
      subproblem.compute_regularized_augmented_matrix_sparsity(view(this->matrix_row_indices), view(this->matrix_column_indices),
         this->solver_indexing);
      this->matrix_values.resize(this->number_nonzeros);
      this->rhs.resize(this->dimension);
      this->solution.resize(this->dimension);
   }

   void COOLinearSystem::initialize(const NewtonSystem& newton_system) {
      this->dimension = newton_system.dimension();
      this->number_nonzeros = newton_system.number_nonzeros();
      this->matrix_row_indices.resize(this->number_nonzeros);
      this->matrix_column_indices.resize(this->number_nonzeros);
      newton_system.compute_sparsity(view(this->matrix_row_indices), view(this->matrix_column_indices), this->solver_indexing);
      this->matrix_values.resize(this->number_nonzeros);
      this->rhs.resize(this->dimension);
      this->solution.resize(this->dimension);
   }

   // the triangular COO representation stores each off-diagonal entry once; duplicate entries are summed
   void COOLinearSystem::compute_matrix_vector_product(const double* vector, double* result, size_t excluded_begin,
         size_t excluded_end) const {
      std::fill(result, result + this->dimension, 0.);
      for (size_t nonzero_index: Range(this->number_nonzeros)) {
         if (excluded_begin <= nonzero_index && nonzero_index < excluded_end) {
            continue;
         }
         const size_t row_index = static_cast<size_t>(this->matrix_row_indices[nonzero_index] - this->solver_indexing);
         const size_t column_index = static_cast<size_t>(this->matrix_column_indices[nonzero_index] - this->solver_indexing);
         const double entry = this->matrix_values[nonzero_index];
         result[row_index] += entry * vector[column_index];
         if (row_index != column_index) {
            result[column_index] += entry * vector[row_index];
         }
      }
   }

   // the quadratic form is not computed (0): the curvature term of the predicted objective reduction is only used by the
   // trust-region method, which does not support the direct linear solvers. The line search uses first-order models
   double COOLinearSystem::compute_hessian_quadratic_form(const Subproblem& /*subproblem*/, const Iterate& /*current_iterate*/,
         const Vector<double>& /*vector*/) const {
      return 0.;
   }
} // namespace