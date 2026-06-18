// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cassert>
#include <numeric>
#include <stdexcept>
#include "HiGHSWorkspace.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "linear_algebra/Indexing.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "symbolic/Range.hpp"

namespace uno {
   void HiGHSWorkspace::allocate_basics(size_t number_variables, size_t number_constraints) {
      this->model.lp_.num_col_ = static_cast<HighsInt>(number_variables);
      this->model.lp_.num_row_ = static_cast<HighsInt>(number_constraints);
      this->model.lp_.sense_ = ObjSense::kMinimize;
      this->model.lp_.offset_ = 0.;

      this->constraints.resize(number_constraints);
      this->linear_objective.resize(number_variables);
      // the linear part of the objective is a dense vector
      this->model.lp_.col_cost_.resize(number_variables);
      // variable bounds
      this->model.lp_.col_lower_.resize(number_variables);
      this->model.lp_.col_upper_.resize(number_variables);
      // constraint bounds
      this->model.lp_.row_lower_.resize(number_constraints);
      this->model.lp_.row_upper_.resize(number_constraints);
   }

   void HiGHSWorkspace::initialize_memory(const Subproblem& subproblem) {
      this->allocate_basics(subproblem.number_variables, subproblem.number_constraints);
      // Jacobian sparsity
      this->compute_jacobian_sparsity(subproblem);
      // Lagrangian Hessian sparsity
      this->compute_hessian_sparsity(subproblem);
   }

   void HiGHSWorkspace::set_from_coo(size_t number_variables, size_t number_constraints, const Vector<double>& linear_objective,
         const Vector<uno_int>& jacobian_row_indices, const Vector<uno_int>& jacobian_column_indices,
         const Vector<double>& jacobian_values,
         const Vector<uno_int>& hessian_row_indices, const Vector<uno_int>& hessian_column_indices,
         const Vector<double>& hessian_values) {
      const size_t number_jacobian_nonzeros = jacobian_values.size();
      const size_t number_hessian_nonzeros = hessian_values.size();
      this->allocate_basics(number_variables, number_constraints);

      // dense objective gradient
      for (size_t variable_index: Range(number_variables)) {
         this->model.lp_.col_cost_[variable_index] = linear_objective[variable_index];
      }

      // constraint Jacobian: copy the COO sparsity + values, build the CSC sparsity, scatter the values
      this->jacobian_row_indices.resize(number_jacobian_nonzeros);
      this->jacobian_column_indices.resize(number_jacobian_nonzeros);
      this->jacobian_values.resize(number_jacobian_nonzeros);
      for (size_t nonzero_index: Range(number_jacobian_nonzeros)) {
         this->jacobian_row_indices[nonzero_index] = jacobian_row_indices[nonzero_index];
         this->jacobian_column_indices[nonzero_index] = jacobian_column_indices[nonzero_index];
         this->jacobian_values[nonzero_index] = jacobian_values[nonzero_index];
      }
      this->build_csc_jacobian_from_coo(number_variables, number_constraints);
      this->scatter_jacobian_values();

      // Lagrangian Hessian: same conversion, or no Hessian (LP) if empty
      if (0 < number_hessian_nonzeros) {
         this->hessian_row_indices.resize(number_hessian_nonzeros);
         this->hessian_column_indices.resize(number_hessian_nonzeros);
         this->hessian_values.resize(number_hessian_nonzeros);
         for (size_t nonzero_index: Range(number_hessian_nonzeros)) {
            this->hessian_row_indices[nonzero_index] = hessian_row_indices[nonzero_index];
            this->hessian_column_indices[nonzero_index] = hessian_column_indices[nonzero_index];
            this->hessian_values[nonzero_index] = hessian_values[nonzero_index];
         }
         this->build_csc_hessian_from_coo(number_variables);
         this->scatter_hessian_values();
      }
      else {
         this->model.hessian_.dim_ = 0;
      }
   }

   double HiGHSWorkspace::compute_hessian_quadratic_form(const Subproblem& /*subproblem*/, const Vector<double>& vector) const {
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

   void HiGHSWorkspace::evaluate_functions(Statistics& statistics, const Subproblem& subproblem,
         Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) {
      // evaluate the functions based on warmstart information
      if (warmstart_information.new_iterate) {
         for (size_t index: Range(subproblem.number_variables)) {
            this->model.lp_.col_cost_[index] = 0.;
         }
         subproblem.problem.evaluate_objective_gradient(subproblem.current_iterate, this->model.lp_.col_cost_.data(),
            current_evaluations);
         subproblem.problem.evaluate_constraints(subproblem.current_iterate, this->constraints.data(), current_evaluations);
         this->evaluate_jacobian(subproblem.problem, subproblem.current_iterate.primals, current_evaluations);
         // evaluate the Hessian and regularize it
         subproblem.evaluate_lagrangian_hessian(statistics, this->hessian_values.data());
         // copy the Hessian with permutation into this->model.hessian_.value_
         this->scatter_hessian_values();
         subproblem.regularize_lagrangian_hessian(statistics, this->model.hessian_.value_.data());
      }
   }

   void HiGHSWorkspace::compute_jacobian_sparsity(const Subproblem& subproblem) {
      const size_t number_jacobian_nonzeros = subproblem.number_jacobian_nonzeros();
      // compute the COO sparsity pattern (column-wise)
      this->jacobian_row_indices.resize(number_jacobian_nonzeros);
      this->jacobian_column_indices.resize(number_jacobian_nonzeros);
      this->jacobian_values.resize(number_jacobian_nonzeros);
      subproblem.compute_jacobian_sparsity(this->jacobian_row_indices.data(), this->jacobian_column_indices.data(), 0, 0,
         Indexing::C_indexing, MatrixOrder::COLUMN_MAJOR);
      // convert COO -> HiGHS' CSC layout
      this->build_csc_jacobian_from_coo(subproblem.number_variables, subproblem.number_constraints);
   }

   // build HiGHS' CSC Jacobian (a_matrix_) and the sorting permutation from the COO arrays already present in
   // jacobian_row_indices/jacobian_column_indices. Shared by the Subproblem path and the data-driven path.
   void HiGHSWorkspace::build_csc_jacobian_from_coo(size_t number_variables, size_t number_constraints) {
      (void) number_constraints;
      const size_t number_jacobian_nonzeros = this->jacobian_row_indices.size();
      // column-wise constraint Jacobian
      this->model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
      this->model.lp_.a_matrix_.index_.resize(number_jacobian_nonzeros); // constraint indices
      this->model.lp_.a_matrix_.start_.resize(number_variables + 1);
      this->model.lp_.a_matrix_.value_.resize(number_jacobian_nonzeros);

      // HiGHS requires a CSC Jacobian: the entries should be in increasing column indices.
      this->jacobian_permutation_vector.resize(number_jacobian_nonzeros);
      std::iota(this->jacobian_permutation_vector.begin(), this->jacobian_permutation_vector.end(), 0);
      std::sort(this->jacobian_permutation_vector.begin(), this->jacobian_permutation_vector.end(),
          [&](const size_t& i, const size_t& j) {
             if (this->jacobian_column_indices[i] < this->jacobian_column_indices[j]) {
                return true;
             }
             // within a given column, have the row indices in increasing order
             else if (this->jacobian_column_indices[i] == this->jacobian_column_indices[j]) {
               return (this->jacobian_row_indices[i] < this->jacobian_row_indices[j]);
             }
             return false;
          }
      );

      this->model.lp_.a_matrix_.start_[0] = 0;
      int current_variable = 0;
      for (size_t jacobian_nonzero_index: Range(number_jacobian_nonzeros)) {
         const size_t permuted_nonzero_index = this->jacobian_permutation_vector[jacobian_nonzero_index];
         // constraint index is used as is
         const HighsInt constraint_index = static_cast<HighsInt>(this->jacobian_row_indices[permuted_nonzero_index]);
         this->model.lp_.a_matrix_.index_[jacobian_nonzero_index] = constraint_index;

         // variable index is used to build the pointers to the column starts
         const uno_int variable_index = this->jacobian_column_indices[permuted_nonzero_index];
         assert(current_variable <= variable_index);
         while (current_variable < variable_index) {
            ++current_variable;
            this->model.lp_.a_matrix_.start_[static_cast<size_t>(current_variable)] = static_cast<HighsInt>(jacobian_nonzero_index);
         }
      }
      // fill the remaining (trailing) empty columns
      while (current_variable < static_cast<int>(number_variables)) {
         ++current_variable;
         this->model.lp_.a_matrix_.start_[static_cast<size_t>(current_variable)] = static_cast<HighsInt>(number_jacobian_nonzeros);
      }
   }

   // column-wise lower triangular Lagrangian Hessian
   void HiGHSWorkspace::compute_hessian_sparsity(const Subproblem& subproblem) {
      if (!subproblem.has_hessian_matrix()) {
         throw std::runtime_error("The subproblem does not have an explicit Hessian matrix and cannot be solved with HiGHS");
      }
      const size_t number_regularized_hessian_nonzeros = subproblem.number_regularized_hessian_nonzeros();
      // get the Hessian sparsity in COO format
      this->hessian_row_indices.resize(number_regularized_hessian_nonzeros);
      this->hessian_column_indices.resize(number_regularized_hessian_nonzeros);
      this->hessian_values.resize(number_regularized_hessian_nonzeros);
      subproblem.compute_regularized_hessian_sparsity(this->hessian_row_indices.data(),
         this->hessian_column_indices.data(), Indexing::C_indexing);
      // convert COO -> HiGHS' lower-triangular CSC layout
      this->build_csc_hessian_from_coo(subproblem.number_variables);
   }

   // build HiGHS' lower-triangular CSC Hessian (model.hessian_) and the sorting permutation from the COO arrays
   // already present in hessian_row_indices/hessian_column_indices. Shared by both build paths.
   void HiGHSWorkspace::build_csc_hessian_from_coo(size_t number_variables) {
      const size_t number_hessian_nonzeros = this->hessian_row_indices.size();
      this->model.hessian_.dim_ = static_cast<HighsInt>(number_variables);
      this->model.hessian_.format_ = HessianFormat::kTriangular;
      this->model.hessian_.index_.resize(number_hessian_nonzeros);
      this->model.hessian_.start_.resize(number_variables + 1);
      this->model.hessian_.value_.resize(number_hessian_nonzeros);

      // HiGHS requires a lower-triangular CSC Hessian: the entries should be in increasing column indices.
      this->hessian_permutation_vector.resize(number_hessian_nonzeros);
      std::iota(this->hessian_permutation_vector.begin(), this->hessian_permutation_vector.end(), 0);
      std::sort(this->hessian_permutation_vector.begin(), this->hessian_permutation_vector.end(),
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
      for (size_t hessian_nonzero_index: Range(number_hessian_nonzeros)) {
         const size_t permuted_nonzero_index = this->hessian_permutation_vector[hessian_nonzero_index];
         // row index
         const HighsInt row_index = static_cast<HighsInt>(this->hessian_row_indices[permuted_nonzero_index]);
         this->model.hessian_.index_[hessian_nonzero_index] = row_index;

         // column index
         const uno_int column_index = this->hessian_column_indices[permuted_nonzero_index];
         assert(current_column <= column_index);
         while (current_column < column_index) {
            ++current_column;
            this->model.hessian_.start_[static_cast<size_t>(current_column)] = static_cast<HighsInt>(hessian_nonzero_index);
         }
      }
      // fill the remaining empty columns
      while (current_column < static_cast<int>(number_variables)) {
         ++current_column;
         this->model.hessian_.start_[static_cast<size_t>(current_column)] = static_cast<HighsInt>(number_hessian_nonzeros);
      }
   }

   void HiGHSWorkspace::scatter_jacobian_values() {
      // copy the Jacobian values with permutation into this->model.lp_.a_matrix_.value_
      const size_t number_jacobian_nonzeros = this->jacobian_values.size();
      for (size_t nonzero_index: Range(number_jacobian_nonzeros)) {
         const size_t permuted_nonzero_index = this->jacobian_permutation_vector[nonzero_index];
         this->model.lp_.a_matrix_.value_[nonzero_index] = this->jacobian_values[permuted_nonzero_index];
      }
   }

   void HiGHSWorkspace::scatter_hessian_values() {
      // copy the Hessian values with permutation into this->model.hessian_.value_
      const size_t number_hessian_nonzeros = this->hessian_values.size();
      for (size_t nonzero_index: Range(number_hessian_nonzeros)) {
         const size_t permuted_nonzero_index = this->hessian_permutation_vector[nonzero_index];
         this->model.hessian_.value_[nonzero_index] = this->hessian_values[permuted_nonzero_index];
      }
   }

   void HiGHSWorkspace::evaluate_jacobian(const OptimizationProblem& problem, const Vector<double>& primals,
         Evaluations& evaluations) {
      problem.evaluate_jacobian(primals, this->jacobian_values.data(), evaluations);
      this->scatter_jacobian_values();
   }
} // namespace
