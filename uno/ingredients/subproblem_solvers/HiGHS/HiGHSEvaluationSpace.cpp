// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include "HiGHSEvaluationSpace.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/Indexing.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/WarmstartInformation.hpp"

namespace uno {
   void HiGHSEvaluationSpace::initialize_memory(const Subproblem& subproblem) {
      this->model.lp_.num_col_ = static_cast<HighsInt>(subproblem.number_variables);
      this->model.lp_.num_row_ = static_cast<HighsInt>(subproblem.number_constraints);

      // determine whether the subproblem has curvature. For the moment, HiGHS can only solve LPs
      this->constraints.resize(subproblem.number_constraints);
      this->linear_objective.resize(subproblem.number_variables);
      this->model.lp_.sense_ = ObjSense::kMinimize;
      this->model.lp_.offset_ = 0.;
      // the linear part of the objective is a dense vector
      this->model.lp_.col_cost_.resize(subproblem.number_variables);
      // variable bounds
      this->model.lp_.col_lower_.resize(subproblem.number_variables);
      this->model.lp_.col_upper_.resize(subproblem.number_variables);
      // constraint bounds
      this->model.lp_.row_lower_.resize(subproblem.number_constraints);
      this->model.lp_.row_upper_.resize(subproblem.number_constraints);

      // column-wise constraint Jacobian
      this->model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
      const size_t number_jacobian_nonzeros = subproblem.number_jacobian_nonzeros();
      // compute the COO sparsity pattern
      this->jacobian_row_indices.resize(number_jacobian_nonzeros);
      this->jacobian_column_indices.resize(number_jacobian_nonzeros);
      subproblem.compute_constraint_jacobian_sparsity(this->jacobian_row_indices.data(),
         this->jacobian_column_indices.data(), Indexing::C_indexing, MatrixOrder::COLUMN_MAJOR);
      // HiGHS matrix in CSC format (variable after variable)
      this->model.lp_.a_matrix_.index_.resize(number_jacobian_nonzeros); // constraint indices
      this->model.lp_.a_matrix_.start_.resize(subproblem.number_variables + 1);
      this->model.lp_.a_matrix_.value_.resize(number_jacobian_nonzeros);
      int current_variable = 0;
      for (size_t jacobian_nonzero_index: Range(number_jacobian_nonzeros)) {
         // constraint index is used as is
         const HighsInt constraint_index = static_cast<HighsInt>(this->jacobian_row_indices[jacobian_nonzero_index]);
         this->model.lp_.a_matrix_.index_[jacobian_nonzero_index] = constraint_index;

         // variable index is used to build the pointers to the column starts
         const int variable_index = this->jacobian_column_indices[jacobian_nonzero_index];
         assert(current_variable <= variable_index);
         while (current_variable < variable_index) {
            ++current_variable;
            this->model.lp_.a_matrix_.start_[static_cast<size_t>(current_variable)] = static_cast<HighsInt>(jacobian_nonzero_index);
         }
      }
      this->model.lp_.a_matrix_.start_[subproblem.number_variables] = static_cast<HighsInt>(number_jacobian_nonzeros);

      // Lagrangian Hessian
      this->compute_hessian_sparsity(subproblem);
   }

   void HiGHSEvaluationSpace::evaluate_constraint_jacobian(const OptimizationProblem& problem, Iterate& iterate) {
      problem.evaluate_constraint_jacobian(iterate, this->model.lp_.a_matrix_.value_.data());
   }

   void HiGHSEvaluationSpace::compute_constraint_jacobian_vector_product(const Vector<double>& vector, Vector<double>& result) const {
      result.fill(0.);
      const size_t number_constraint_jacobian_nonzeros = this->jacobian_row_indices.size();
      for (size_t nonzero_index: Range(number_constraint_jacobian_nonzeros)) {
         const size_t constraint_index = static_cast<size_t>(this->jacobian_row_indices[nonzero_index]);
         const size_t variable_index = static_cast<size_t>(this->jacobian_column_indices[nonzero_index]);
         const double derivative = this->model.lp_.a_matrix_.value_[nonzero_index];

         // a safeguard to make sure we take only the correct part of the Jacobian
         if (variable_index < vector.size() && constraint_index < result.size()) {
            result[constraint_index] += derivative * vector[variable_index];
         }
      }
   }

   void HiGHSEvaluationSpace::compute_constraint_jacobian_transposed_vector_product(const Vector<double>& vector, Vector<double>& result) const {
      result.fill(0.);
      const size_t number_constraint_jacobian_nonzeros = this->jacobian_row_indices.size();
      for (size_t nonzero_index: Range(number_constraint_jacobian_nonzeros)) {
         const size_t constraint_index = static_cast<size_t>(this->jacobian_row_indices[nonzero_index]);
         const size_t variable_index = static_cast<size_t>(this->jacobian_column_indices[nonzero_index]);
         const double derivative = this->model.lp_.a_matrix_.value_[nonzero_index];
         assert(constraint_index < vector.size());
         assert(variable_index < result.size());

         result[variable_index] += derivative * vector[constraint_index];
      }
   }

   double HiGHSEvaluationSpace::compute_hessian_quadratic_product(const Vector<double>& vector) const {
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

   void HiGHSEvaluationSpace::evaluate_functions(Statistics& statistics, const Subproblem& subproblem,
         const WarmstartInformation& warmstart_information) {
      // evaluate the functions based on warmstart information
      if (warmstart_information.objective_changed) {
         subproblem.problem.evaluate_objective_gradient(subproblem.current_iterate, this->model.lp_.col_cost_.data());
      }
      if (warmstart_information.constraints_changed) {
         subproblem.problem.evaluate_constraints(subproblem.current_iterate, this->constraints);
         this->evaluate_constraint_jacobian(subproblem.problem, subproblem.current_iterate);
      }
      // evaluate the Hessian and regularize it
      if (warmstart_information.objective_changed || warmstart_information.constraints_changed) {
         subproblem.evaluate_lagrangian_hessian(statistics, this->hessian_values.data());
         // copy the Hessian with permutation into this->model.hessian_.value_
         for (size_t nonzero_index: Range(subproblem.number_regularized_hessian_nonzeros())) {
            const size_t permutated_nonzero_index = this->permutation_vector[nonzero_index];
            this->model.hessian_.value_[nonzero_index] = this->hessian_values[permutated_nonzero_index];
         }
         subproblem.regularize_lagrangian_hessian(statistics, this->model.hessian_.value_.data());
      }
   }

   // column-wise lower triangular Lagrangian Hessian
   void HiGHSEvaluationSpace::compute_hessian_sparsity(const Subproblem& subproblem) {
      if (!subproblem.has_hessian_matrix()) {
         throw std::runtime_error("The subproblem does not have an explicit Hessian matrix and cannot be solved with HiGHS");
      }
      const size_t number_regularized_hessian_nonzeros = subproblem.number_regularized_hessian_nonzeros();
      this->model.hessian_.dim_ = static_cast<HighsInt>(subproblem.number_variables);
      this->model.hessian_.format_ = HessianFormat::kTriangular;
      this->model.hessian_.index_.resize(number_regularized_hessian_nonzeros);
      this->model.hessian_.start_.resize(subproblem.number_variables + 1);
      this->model.hessian_.value_.resize(number_regularized_hessian_nonzeros);

      // get the Jacobian sparsity in COO format
      this->hessian_row_indices.resize(number_regularized_hessian_nonzeros);
      this->hessian_column_indices.resize(number_regularized_hessian_nonzeros);
      subproblem.compute_regularized_hessian_sparsity(this->hessian_row_indices.data(),
         this->hessian_column_indices.data(), Indexing::C_indexing);

      // HiGHS requires a lower-triangular CSC Hessian: the entries should be in increasing column indices.
      // Since the COO format does not require this, we need to convert from COO to CSC by permutating the entries. To
      // this end, we compute a permutation vector once and for all that contains the correct ordering of terms.
      // The permutation vector is initially filled with [0, 1, ..., number_regularized_hessian_nonzeros-1]
      this->permutation_vector.resize(number_regularized_hessian_nonzeros);
      std::iota(this->permutation_vector.begin(), this->permutation_vector.end(), 0);
      // sort the permutation vector such that the column indices are in increasing order
      // see https://stackoverflow.com/questions/17554242/how-to-obtain-the-index-permutation-after-the-sorting
      std::sort(this->permutation_vector.begin(), this->permutation_vector.end(),
          [&](const size_t& i, const size_t& j) {
             if (this->hessian_column_indices[i] < this->hessian_column_indices[j]) {
                return true;
             }
             // within a given column, have the row indices in increasing order
             else if (this->hessian_column_indices[i] == this->hessian_column_indices[j]) {
               return (this->hessian_row_indices[i] < this->hessian_row_indices[j]);
             }
             return false;
          }
      );

      // copy the COO format into HiGHS' CSC format
      this->model.hessian_.start_[0] = 0;
      int current_column = 0;
      for (size_t hessian_nonzero_index: Range(number_regularized_hessian_nonzeros)) {
         const size_t permutated_nonzero_index = this->permutation_vector[hessian_nonzero_index];
         // row index
         const HighsInt row_index = static_cast<HighsInt>(this->hessian_row_indices[permutated_nonzero_index]);
         this->model.hessian_.index_[hessian_nonzero_index] = row_index;

         // column index
         const int column_index = this->hessian_column_indices[permutated_nonzero_index];
         assert(current_column <= column_index);
         while (current_column < column_index) {
            ++current_column;
            this->model.hessian_.start_[static_cast<size_t>(current_column)] = static_cast<HighsInt>(hessian_nonzero_index);
         }
      }
      // fill the remaining empty columns
      while (current_column < static_cast<int>(subproblem.number_variables)) {
         ++current_column;
         this->model.hessian_.start_[static_cast<size_t>(current_column)] = static_cast<HighsInt>(number_regularized_hessian_nonzeros);
      }

      // the Hessian will be evaluated in this vector, and copied with permutation into this->model.hessian_.value_
      this->hessian_values.resize(number_regularized_hessian_nonzeros);
   }
} // namespace