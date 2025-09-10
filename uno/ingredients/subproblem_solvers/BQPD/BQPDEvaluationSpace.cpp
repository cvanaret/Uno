// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <numeric>
#include "BQPDEvaluationSpace.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/Indexing.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/WarmstartInformation.hpp"

namespace uno {
   void BQPDEvaluationSpace::initialize(const Subproblem& subproblem) {
      this->constraints.resize(subproblem.number_constraints);

      // Jacobian + objective gradient
      this->gradients.resize(subproblem.number_variables + subproblem.number_jacobian_nonzeros());
      this->gradient_sparsity.resize(1 + subproblem.number_variables + subproblem.number_jacobian_nonzeros() +
         1 + subproblem.number_constraints + 1);
      // save sparsity patterns of objective gradient and constraint Jacobian into BQPD workspace
      this->compute_gradients_sparsity(subproblem);

      // allocate an explicit Hessian matrix if:
      // - the Hessian is not positive definite and must be regularized, or
      // - the Hessian model only has an explicit representation
      if ((!subproblem.is_hessian_positive_definite() && subproblem.performs_primal_regularization()) ||
            !subproblem.has_hessian_operator()) {
         const size_t number_regularized_hessian_nonzeros = subproblem.number_regularized_hessian_nonzeros();
         this->hessian_row_indices.resize(number_regularized_hessian_nonzeros);
         this->hessian_column_indices.resize(number_regularized_hessian_nonzeros);
         this->hessian_values.resize(number_regularized_hessian_nonzeros);
         subproblem.compute_regularized_hessian_sparsity(this->hessian_row_indices.data(),
            this->hessian_column_indices.data(), Indexing::C_indexing);
      }
   }

   void BQPDEvaluationSpace::evaluate_constraint_jacobian(const OptimizationProblem& problem, Iterate& iterate) {
      problem.evaluate_constraint_jacobian(iterate, this->jacobian_values.data());

      // copy the Jacobian with permutation into &this->gradients[subproblem.number_variables]
      for (size_t nonzero_index: Range(problem.number_jacobian_nonzeros())) {
         const size_t permutated_nonzero_index = this->permutation_vector[nonzero_index];
         this->gradients[problem.number_variables + nonzero_index] = this->jacobian_values[permutated_nonzero_index];
      }
   }

   void BQPDEvaluationSpace::compute_constraint_jacobian_vector_product(const Vector<double>& vector, Vector<double>& result) const {
      result.fill(0.);
      const size_t number_constraint_jacobian_nonzeros = this->jacobian_row_indices.size();
      for (size_t nonzero_index: Range(number_constraint_jacobian_nonzeros)) {
         const size_t constraint_index = static_cast<size_t>(this->jacobian_row_indices[nonzero_index]);
         const size_t variable_index = static_cast<size_t>(this->jacobian_column_indices[nonzero_index]);
         const double derivative = this->jacobian_values[nonzero_index];

         // a safeguard to make sure we take only the correct part of the Jacobian
         if (variable_index < vector.size() && constraint_index < result.size()) {
            result[constraint_index] += derivative * vector[variable_index];
         }
      }
   }

   void BQPDEvaluationSpace::compute_constraint_jacobian_transposed_vector_product(const Vector<double>& vector, Vector<double>& result) const {
      result.fill(0.);
      const size_t number_constraint_jacobian_nonzeros = this->jacobian_row_indices.size();
      for (size_t nonzero_index: Range(number_constraint_jacobian_nonzeros)) {
         const size_t constraint_index = static_cast<size_t>(this->jacobian_row_indices[nonzero_index]);
         const size_t variable_index = static_cast<size_t>(this->jacobian_column_indices[nonzero_index]);
         const double derivative = this->jacobian_values[nonzero_index];
         assert(constraint_index < vector.size());
         assert(variable_index < result.size());

         result[variable_index] += derivative * vector[constraint_index];
      }
   }

   double BQPDEvaluationSpace::compute_hessian_quadratic_product(const Vector<double>& vector) const {
      double quadratic_product = 0.;
      const size_t number_hessian_nonzeros = this->hessian_values.size();
      for (size_t nonzero_index: Range(number_hessian_nonzeros)) {
         const size_t row_index = static_cast<size_t>(this->hessian_row_indices[nonzero_index]);
         const size_t column_index = static_cast<size_t>(this->hessian_column_indices[nonzero_index]);
         const double entry = this->hessian_values[nonzero_index];
         assert(row_index < vector.size());
         assert(column_index < vector.size());

         const double factor = (row_index != column_index) ? 2. : 1.;
         quadratic_product += factor * entry * vector[row_index] * vector[column_index];
      }
      return quadratic_product;
   }

   void BQPDEvaluationSpace::evaluate_functions(const OptimizationProblem& problem, Iterate& current_iterate,
         const WarmstartInformation& warmstart_information) {
      // evaluate the functions based on warmstart information
      // gradients is a concatenation of the dense objective gradient and the sparse Jacobian
      if (warmstart_information.objective_changed) {
         problem.evaluate_objective_gradient(current_iterate, this->gradients.data());
      }
      if (warmstart_information.constraints_changed) {
         problem.evaluate_constraints(current_iterate, this->constraints);
         this->evaluate_constraint_jacobian(problem, current_iterate);
      }
      if (warmstart_information.objective_changed || warmstart_information.constraints_changed) {
         this->evaluate_hessian = true;
      }
   }

   // objective gradient + row-major constraint Jacobian
   void BQPDEvaluationSpace::compute_gradients_sparsity(const Subproblem& subproblem) {
      const size_t number_jacobian_nonzeros = subproblem.number_jacobian_nonzeros();

      // header
      const size_t position_of_row_starts = 1 + subproblem.number_variables + number_jacobian_nonzeros;
      this->gradient_sparsity[0] = static_cast<int>(position_of_row_starts);

      // dense objective gradient
      for (size_t variable_index: Range(subproblem.number_variables)) {
         this->gradient_sparsity[1 + variable_index] = static_cast<int>(variable_index + Indexing::Fortran_indexing);
      }

      this->gradient_sparsity[position_of_row_starts] = 1; // always starts at 1
      this->gradient_sparsity[position_of_row_starts + 1] = static_cast<int>(1 + subproblem.number_variables); // dense objective gradient

      // get the Jacobian sparsity in COO format
      this->jacobian_row_indices.resize(number_jacobian_nonzeros);
      this->jacobian_column_indices.resize(number_jacobian_nonzeros);
      subproblem.compute_constraint_jacobian_sparsity(this->jacobian_row_indices.data(),
         this->jacobian_column_indices.data(), Indexing::C_indexing, MatrixOrder::ROW_MAJOR);

      // BQPD (sparse) requires a (weak) CSR Jacobian: the entries should be in increasing constraint indices.
      // Since the COO format does not require this, we need to convert from COO to CSR by permutating the entries. To
      // this end, we compute a permutation vector once and for all that contains the correct ordering of terms.
      // The permutation vector is initially filled with [0, 1, ..., number_jacobian_nonzeros-1]
      this->permutation_vector.resize(number_jacobian_nonzeros);
      std::iota(this->permutation_vector.begin(), this->permutation_vector.end(), 0);
      // sort the permutation vector such that the row indices (constraints) of the Jacobian sparsity are in increasing order
      // see https://stackoverflow.com/questions/17554242/how-to-obtain-the-index-permutation-after-the-sorting
      std::sort(this->permutation_vector.begin(), this->permutation_vector.end(),
          [&](const size_t& i, const size_t& j) {
             return (this->jacobian_row_indices[i] < this->jacobian_row_indices[j]);
          }
      );

      // copy the COO format into BQPD's CSR format
      int current_constraint = 0;
      for (size_t jacobian_nonzero_index: Range(number_jacobian_nonzeros)) {
         const size_t permutated_nonzero_index = this->permutation_vector[jacobian_nonzero_index];
         // variable index
         const int variable_index = this->jacobian_column_indices[permutated_nonzero_index];
         this->gradient_sparsity[1 + subproblem.number_variables + jacobian_nonzero_index] = variable_index +
            Indexing::Fortran_indexing;

         // constraint index
         const int constraint_index = this->jacobian_row_indices[permutated_nonzero_index];
         assert(current_constraint <= constraint_index);
         while (current_constraint < constraint_index) {
            ++current_constraint;
            this->gradient_sparsity[1 + subproblem.number_variables + number_jacobian_nonzeros + 1 +
               static_cast<size_t>(current_constraint)] = static_cast<int>(subproblem.number_variables +
                  jacobian_nonzero_index) + Indexing::Fortran_indexing;
         }
      }
      // since there cannot be empty rows, we don't need to loop over empty rows like we do for the HiGHS Hessian
      this->gradient_sparsity[1 + subproblem.number_variables + number_jacobian_nonzeros + 1 + subproblem.number_constraints] =
         static_cast<int>(subproblem.number_variables + number_jacobian_nonzeros) + Indexing::Fortran_indexing;

      // the Jacobian will be evaluated in this vector, and copied with permutation into this->gradients
      this->jacobian_values.resize(number_jacobian_nonzeros);
   }
} // namespace
