// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "BQPDEvaluationSpace.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/WarmstartInformation.hpp"

namespace uno {
   void BQPDEvaluationSpace::evaluate_constraint_jacobian(const Subproblem& subproblem) {
      subproblem.evaluate_constraint_jacobian(this->jacobian_values.data());

      // copy the Jacobian with permutation into &this->gradients[subproblem.number_variables]
      for (size_t nonzero_index: Range(subproblem.number_jacobian_nonzeros())) {
         const size_t permutated_nonzero_index = this->permutation_vector[nonzero_index];
         this->gradients[subproblem.number_variables + nonzero_index] = this->jacobian_values[permutated_nonzero_index];
      }
   }

   void BQPDEvaluationSpace::compute_constraint_jacobian_vector_product(const Vector<double>& vector, Vector<double>& result) const {
      result.fill(0.);
      const size_t number_constraint_jacobian_nonzeros = this->jacobian_row_indices.size();
      for (size_t nonzero_index: Range(number_constraint_jacobian_nonzeros)) {
         const size_t constraint_index = this->jacobian_row_indices[nonzero_index];
         const size_t variable_index = this->jacobian_column_indices[nonzero_index];
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
         const size_t constraint_index = this->jacobian_row_indices[nonzero_index];
         const size_t variable_index = this->jacobian_column_indices[nonzero_index];
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
         const size_t row_index = this->hessian_row_indices[nonzero_index];
         const size_t column_index = this->hessian_column_indices[nonzero_index];
         const double entry = this->hessian_values[nonzero_index];
         assert(row_index < vector.size());
         assert(column_index < vector.size());

         const double factor = (row_index != column_index) ? 2. : 1.;
         quadratic_product += factor * entry * vector[row_index] * vector[column_index];
      }
      return quadratic_product;
   }

   void BQPDEvaluationSpace::evaluate_functions(const Subproblem& subproblem, const WarmstartInformation& warmstart_information) {
      // evaluate the functions based on warmstart information
      // gradients is a concatenation of the dense objective gradient and the sparse Jacobian
      if (warmstart_information.objective_changed) {
         subproblem.evaluate_objective_gradient(this->gradients.data());
      }
      if (warmstart_information.constraints_changed) {
         subproblem.evaluate_constraints(this->constraints);
         this->evaluate_constraint_jacobian(subproblem);
      }
      if (warmstart_information.objective_changed || warmstart_information.constraints_changed) {
         this->evaluate_hessian = true;
      }
   }
} // namespace