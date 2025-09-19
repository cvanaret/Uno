// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include "COOEvaluationSpace.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "linear_algebra/COOMatrix.hpp"
#include "linear_algebra/Indexing.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/WarmstartInformation.hpp"

namespace uno {
   void COOEvaluationSpace::initialize_hessian(const Subproblem& subproblem) {
      if (!subproblem.has_hessian_matrix()) {
         throw std::runtime_error("The subproblem does not have an explicit Hessian matrix and cannot be solved with a direct linear solver");
      }
      const size_t dimension = subproblem.number_variables;

      // Hessian
      this->number_hessian_nonzeros = subproblem.number_hessian_nonzeros();
      this->number_matrix_nonzeros = subproblem.number_regularized_hessian_nonzeros();
      this->matrix_row_indices.resize(this->number_matrix_nonzeros);
      this->matrix_column_indices.resize(this->number_matrix_nonzeros);
      // compute the COO sparse representation
      subproblem.compute_regularized_hessian_sparsity(this->matrix_row_indices.data(), this->matrix_column_indices.data(),
         Indexing::Fortran_indexing);
      this->matrix_values.resize(this->number_matrix_nonzeros);
      this->rhs.resize(dimension);
      this->solution.resize(dimension);
   }

   void COOEvaluationSpace::initialize_augmented_system(const Subproblem& subproblem) {
      if (!subproblem.has_hessian_matrix()) {
         throw std::runtime_error("The subproblem does not have an explicit Hessian matrix and cannot be solved with a direct linear solver");
      }
      const size_t dimension = subproblem.number_variables + subproblem.number_constraints;

      // evaluations
      this->objective_gradient.resize(subproblem.number_variables);
      this->constraints.resize(subproblem.number_constraints);

      // Jacobian
      this->number_jacobian_nonzeros = subproblem.number_jacobian_nonzeros();
      this->jacobian_row_indices.resize(this->number_jacobian_nonzeros);
      this->jacobian_column_indices.resize(this->number_jacobian_nonzeros);
      subproblem.compute_constraint_jacobian_sparsity(this->jacobian_row_indices.data(), this->jacobian_column_indices.data(),
         Indexing::C_indexing, MatrixOrder::COLUMN_MAJOR);

      // augmented system
      this->number_hessian_nonzeros = subproblem.number_hessian_nonzeros();
      this->number_matrix_nonzeros = subproblem.number_regularized_augmented_system_nonzeros();
      this->matrix_row_indices.resize(this->number_matrix_nonzeros);
      this->matrix_column_indices.resize(this->number_matrix_nonzeros);
      // compute the COO sparse representation
      subproblem.compute_regularized_augmented_matrix_sparsity(this->matrix_row_indices.data(), this->matrix_column_indices.data(),
         this->jacobian_row_indices.data(), this->jacobian_column_indices.data(), Indexing::Fortran_indexing);
      this->matrix_values.resize(this->number_matrix_nonzeros);
      this->rhs.resize(dimension);
      this->solution.resize(dimension);
   }

   void COOEvaluationSpace::evaluate_constraint_jacobian(const OptimizationProblem& problem, Iterate& iterate) {
      problem.evaluate_constraint_jacobian(iterate, this->matrix_values.data() + this->number_hessian_nonzeros);
   }

   void COOEvaluationSpace::compute_constraint_jacobian_vector_product(const Vector<double>& vector, Vector<double>& result) const {
      result.fill(0.);
      const size_t offset = this->number_hessian_nonzeros;
      for (size_t nonzero_index: Range(this->number_jacobian_nonzeros)) {
         const size_t constraint_index = static_cast<size_t>(this->jacobian_row_indices[nonzero_index]);
         const size_t variable_index = static_cast<size_t>(this->jacobian_column_indices[nonzero_index]);
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
         const size_t constraint_index = static_cast<size_t>(this->jacobian_row_indices[nonzero_index]);
         const size_t variable_index = static_cast<size_t>(this->jacobian_column_indices[nonzero_index]);
         const double derivative = this->matrix_values[offset + nonzero_index];

         if (variable_index < result.size() && constraint_index < vector.size()) {
            result[variable_index] += derivative * vector[constraint_index];
         }
      }
   }

   double COOEvaluationSpace::compute_hessian_quadratic_product(const Vector<double>& /*vector*/) const {
      return 0.;
   }

   void COOEvaluationSpace::set_up_linear_system(Statistics& statistics, const Subproblem& subproblem,
         DirectSymmetricIndefiniteLinearSolver<double>& linear_solver, const WarmstartInformation& warmstart_information) {
      // evaluate the functions at the current iterate
      if (warmstart_information.objective_changed) {
         subproblem.problem.evaluate_objective_gradient(subproblem.current_iterate, this->objective_gradient.data());
      }
      if (warmstart_information.constraints_changed) {
         subproblem.problem.evaluate_constraints(subproblem.current_iterate, this->constraints);
      }

      if (warmstart_information.objective_changed || warmstart_information.constraints_changed) {
         // perform the symbolic analysis once and for all
         if (!this->analysis_performed) {
            DEBUG << "Performing symbolic analysis of the indefinite system\n";
            linear_solver.do_symbolic_analysis();
            this->analysis_performed = true;
         }
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
} // namespace