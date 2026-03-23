// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "COOLinearSolverSparseRepresentation.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/Indexing.hpp"

namespace uno {
   void COOLinearSolverSparseRepresentation::initialize_hessian(const Subproblem& subproblem) {
      // Hessian
      this->dimension = subproblem.number_variables;
      this->number_hessian_nonzeros = subproblem.number_hessian_nonzeros();
      this->number_nonzeros = subproblem.number_regularized_hessian_nonzeros();
      this->matrix_row_indices.resize(this->number_nonzeros);
      this->matrix_column_indices.resize(this->number_nonzeros);
      // compute the COO sparse representation
      subproblem.compute_regularized_hessian_sparsity(this->matrix_row_indices.data(), this->matrix_column_indices.data(),
         Indexing::Fortran_indexing);
      this->matrix_values.resize(this->number_nonzeros);
      const size_t dimension = subproblem.number_variables;
      this->rhs.resize(dimension);
      this->solution.resize(dimension);
   }

   void COOLinearSolverSparseRepresentation::initialize_augmented_system(const Subproblem& subproblem) {
      // Jacobian
      this->number_jacobian_nonzeros = subproblem.number_jacobian_nonzeros();
      this->jacobian_row_indices.resize(this->number_jacobian_nonzeros);
      this->jacobian_column_indices.resize(this->number_jacobian_nonzeros);
      subproblem.compute_jacobian_sparsity(this->jacobian_row_indices.data(), this->jacobian_column_indices.data(),
         Indexing::C_indexing, MatrixOrder::COLUMN_MAJOR);

      // augmented system
      this->dimension = subproblem.number_variables + subproblem.number_constraints;
      this->number_hessian_nonzeros = subproblem.number_hessian_nonzeros();
      this->number_nonzeros = subproblem.number_regularized_augmented_system_nonzeros();
      this->matrix_row_indices.resize(this->number_nonzeros);
      this->matrix_column_indices.resize(this->number_nonzeros);
      // compute the COO sparse representation
      subproblem.compute_regularized_augmented_matrix_sparsity(this->matrix_row_indices.data(), this->matrix_column_indices.data(),
         this->jacobian_row_indices.data(), this->jacobian_column_indices.data(), Indexing::Fortran_indexing);
      this->matrix_values.resize(this->number_nonzeros);
      const size_t dimension = subproblem.number_variables + subproblem.number_constraints;
      this->rhs.resize(dimension);
      this->solution.resize(dimension);
   }

   double COOLinearSolverSparseRepresentation::compute_hessian_quadratic_product(const Subproblem& /*subproblem*/,
         const Vector<double>& /*vector*/) const {
      return 0.;
   }
} // namespace
