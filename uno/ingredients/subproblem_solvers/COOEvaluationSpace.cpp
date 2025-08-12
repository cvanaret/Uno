// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "COOEvaluationSpace.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/COOMatrix.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/WarmstartInformation.hpp"

namespace uno {
   void COOEvaluationSpace::evaluate_constraint_jacobian(const Subproblem& subproblem) {
      subproblem.evaluate_constraint_jacobian(this->matrix_values.data() + this->number_hessian_nonzeros);
   }

   void COOEvaluationSpace::set_up_linear_system(Statistics& statistics, const Subproblem& subproblem,
         DirectSymmetricIndefiniteLinearSolver<double>& linear_solver, const WarmstartInformation& warmstart_information) {
      // evaluate the functions at the current iterate
      if (warmstart_information.objective_changed) {
         subproblem.evaluate_objective_gradient(this->objective_gradient.data());
      }
      if (warmstart_information.constraints_changed) {
         subproblem.evaluate_constraints(this->constraints);
      }

      if (warmstart_information.objective_changed || warmstart_information.constraints_changed) {
         // assemble the augmented matrix
         subproblem.assemble_augmented_matrix(statistics, this->matrix_values.data());
         // regularize the augmented matrix (this calls the analysis and the factorization)
         subproblem.regularize_augmented_matrix(statistics, this->matrix_values.data(),
            subproblem.dual_regularization_factor(), linear_solver);

         // assemble the RHS
         const COOMatrix jacobian{this->jacobian_row_indices.data(), this->jacobian_column_indices.data(),
            this->matrix_values.data() + this->number_hessian_nonzeros};
         subproblem.assemble_augmented_rhs(this->objective_gradient, this->constraints, jacobian, this->rhs);
      }
   }

   void COOEvaluationSpace::compute_constraint_jacobian_vector_product(const Vector<double>& vector, Vector<double>& result) const {
      result.fill(0.);
      const size_t offset = this->number_hessian_nonzeros;
      for (size_t nonzero_index: Range(this->number_jacobian_nonzeros)) {
         const size_t constraint_index = this->jacobian_row_indices[nonzero_index];
         const size_t variable_index = this->jacobian_column_indices[nonzero_index];
         const double derivative = this->matrix_values[offset + nonzero_index];

         if (constraint_index < result.size() && variable_index < vector.size()) {
            result[constraint_index] += derivative * vector[variable_index];
         }
      }
   }

   void COOEvaluationSpace::compute_constraint_jacobian_transposed_vector_product(const Vector<double>& vector, Vector<double>& result) const {
      result.fill(0.);
      const size_t offset = this->number_hessian_nonzeros;
      for (size_t nonzero_index: Range(this->number_jacobian_nonzeros)) {
         const size_t constraint_index = this->jacobian_row_indices[nonzero_index];
         const size_t variable_index = this->jacobian_column_indices[nonzero_index];
         const double derivative = this->matrix_values[offset + nonzero_index];

         if (variable_index < result.size() && constraint_index < vector.size()) {
            result[variable_index] += derivative * vector[constraint_index];
         }
      }
   }

   double COOEvaluationSpace::compute_hessian_quadratic_product(const Vector<double>& /*vector*/) const {
      throw std::runtime_error("COOEvaluationSpace::compute_hessian_quadratic_product not implemented");
   }
} // namespace