// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "PythonModel.hpp"
#include "symbolic/Concatenation.hpp"
#include "Uno.hpp"

namespace uno {
   PythonModel::PythonModel(const PythonUserModel& user_model):
      Model("Python model", static_cast<size_t>(user_model.number_variables), static_cast<size_t>(user_model.number_constraints),
         static_cast<double>(user_model.optimization_sense)),
      user_model(user_model),
      equality_constraints_collection(this->equality_constraints),
      inequality_constraints_collection(this->inequality_constraints) {
      this->partition_constraints(this->equality_constraints, this->inequality_constraints);
   }

   bool PythonModel::has_implicit_hessian_representation() const {
      return false;
   }

   bool PythonModel::has_explicit_hessian_representation() const {
      return true;
   }

   double PythonModel::evaluate_objective(const Vector<double>& x) const {
      if (this->user_model.objective_function != nullptr) {
         std::cout << "PythonModel::evaluate_objective\n";
         return (*this->user_model.objective_function)(wrap_pointer(&x));
      }
      else {
         std::cout << "PythonModel::evaluate_objective has empty pointer\n";
      }
      return 0.;
   }

   void PythonModel::evaluate_constraints(const Vector<double>& x, std::vector<double>& constraints) const {
      if (this->user_model.constraint_functions != nullptr) {
         this->user_model.constraint_functions(x, constraints);
      }
   }

   void PythonModel::evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const {
      if (this->user_model.objective_gradient != nullptr) {
         std::cout << "PythonModel::evaluate_objective_gradient\n";
         (*this->user_model.objective_gradient)(wrap_pointer(&x), wrap_pointer(&gradient));
      }
   }

   void PythonModel::compute_constraint_jacobian_sparsity(int* row_indices, int* column_indices, int solver_indexing,
         MatrixOrder /*matrix_format*/) const {
      // copy the indices of the user sparsity patterns to the Uno vectors
      std::copy_n(this->user_model.jacobian_row_indices.data(), static_cast<size_t>(this->user_model.number_jacobian_nonzeros), row_indices);
      std::copy_n(this->user_model.jacobian_column_indices.data(), static_cast<size_t>(this->user_model.number_jacobian_nonzeros), column_indices);
      // TODO matrix_format

      // handle the solver indexing
      if (this->user_model.base_indexing != solver_indexing) {
         const int indexing_difference = solver_indexing - this->user_model.base_indexing;
         for (size_t index: Range(static_cast<size_t>(this->user_model.number_jacobian_nonzeros))) {
            row_indices[index] += indexing_difference;
            column_indices[index] += indexing_difference;
         }
      }
   }

   void PythonModel::compute_hessian_sparsity(int* row_indices, int* column_indices, int solver_indexing) const {
      // copy the indices of the user sparsity patterns to the Uno vectors
      std::copy_n(this->user_model.hessian_row_indices.data(), static_cast<size_t>(this->user_model.number_hessian_nonzeros), row_indices);
      std::copy_n(this->user_model.hessian_column_indices.data(), static_cast<size_t>(this->user_model.number_hessian_nonzeros), column_indices);

      // handle the solver indexing
      if (this->user_model.base_indexing != solver_indexing) {
         const int indexing_difference = solver_indexing - this->user_model.base_indexing;
         for (size_t index: Range(static_cast<size_t>(this->user_model.number_hessian_nonzeros))) {
            row_indices[index] += indexing_difference;
            column_indices[index] += indexing_difference;
         }
      }
   }

   void PythonModel::evaluate_constraint_jacobian(const Vector<double>& x, double* jacobian_values) const {
      if (this->user_model.constraint_jacobian != nullptr) {
         this->user_model.constraint_jacobian(x, jacobian_values);
      }
   }

   void PythonModel::evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         double* hessian_values) const {
      if (this->user_model.lagrangian_hessian != nullptr) {
         this->user_model.lagrangian_hessian(x, objective_multiplier, multipliers, hessian_values);
      }
   }

   void PythonModel::compute_hessian_vector_product(const double* /*vector*/, double /*objective_multiplier*/,
         const Vector<double>& /*multipliers*/, double* /*result*/) const {
      throw std::runtime_error("PythonModel::compute_hessian_vector_product not implemented yet");
   }

   double PythonModel::variable_lower_bound(size_t variable_index) const {
      return this->user_model.variables_lower_bounds[variable_index];
   }

   double PythonModel::variable_upper_bound(size_t variable_index) const {
      return this->user_model.variables_upper_bounds[variable_index];
   }

   const SparseVector<size_t>& PythonModel::get_slacks() const {
      return this->slacks;
   }

   const Vector<size_t>& PythonModel::get_fixed_variables() const {
      return this->fixed_variables;
   }

   double PythonModel::constraint_lower_bound(size_t constraint_index) const {
      return this->user_model.constraints_lower_bounds[constraint_index];
   }

   double PythonModel::constraint_upper_bound(size_t constraint_index) const {
      return this->user_model.constraints_upper_bounds[constraint_index];
   }

   const Collection<size_t>& PythonModel::get_equality_constraints() const {
      return this->equality_constraints_collection;
   }

   const Collection<size_t>& PythonModel::get_inequality_constraints() const {
      return this->inequality_constraints_collection;
   }

   const Collection<size_t>& PythonModel::get_linear_constraints() const {
      return this->linear_constraints;
   }

   void PythonModel::initial_primal_point(Vector<double>& x) const {
      if (this->user_model.initial_primal_iterate.has_value()) {
         std::copy_n(this->user_model.initial_primal_iterate->begin(), this->user_model.number_variables, x.begin());
      }
      else {
         x.fill(0.);
      }
   }

   void PythonModel::initial_dual_point(Vector<double>& multipliers) const {
      if (this->user_model.initial_dual_iterate.has_value()) {
         std::copy_n(this->user_model.initial_dual_iterate->begin(), this->user_model.number_constraints, multipliers.begin());
      }
      else {
         multipliers.fill(0.);
      }
   }

   void PythonModel::postprocess_solution(Iterate& /*iterate*/, IterateStatus /*iterate_status*/) const {
      // do nothing
   }

   size_t PythonModel::number_jacobian_nonzeros() const {
      return static_cast<size_t>(this->user_model.number_jacobian_nonzeros);
   }

   size_t PythonModel::number_hessian_nonzeros() const {
      return static_cast<size_t>(this->user_model.number_hessian_nonzeros);
   }
} // namespace