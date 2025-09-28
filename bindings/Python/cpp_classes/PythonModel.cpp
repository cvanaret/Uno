// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "PythonModel.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "symbolic/Concatenation.hpp"
#include "Uno.hpp"

namespace uno {
   PythonModel::PythonModel(const PythonUserModel& user_model):
      Model("Python model", static_cast<size_t>(user_model.number_variables), static_cast<size_t>(user_model.number_constraints),
         static_cast<double>(user_model.optimization_sense)),
         user_model(user_model),
         equality_constraints_collection(this->equality_constraints),
         inequality_constraints_collection(this->inequality_constraints) {
      // find fixed variables
      this->find_fixed_variables(this->fixed_variables);

      // partition equality/inequality constraints
      this->partition_constraints(this->equality_constraints, this->inequality_constraints);
   }

   bool PythonModel::has_jacobian_operator() const {
      return (this->user_model.jacobian_operator != nullptr);
   }

   bool PythonModel::has_jacobian_transposed_operator() const {
      return (this->user_model.jacobian_transposed_operator != nullptr);
   }

   bool PythonModel::has_hessian_operator() const {
      return (this->user_model.lagrangian_hessian_operator != nullptr);
   }

   bool PythonModel::has_hessian_matrix() const {
      return (this->user_model.lagrangian_hessian != nullptr);
   }

   double PythonModel::evaluate_objective(const Vector<double>& x) const {
      double objective_value = 0.;
      if (this->user_model.objective_function.has_value()) {
         const py::object user_data = this->user_model.user_data.has_value() ? *this->user_model.user_data : py::cast(nullptr);
         const int32_t return_code = (*this->user_model.objective_function)(static_cast<int32_t>(this->number_variables),
            x, &objective_value, user_data);
         if (0 < return_code) {
            throw FunctionEvaluationError();
         }
         objective_value *= this->optimization_sense;
      }
      return objective_value;
   }

   void PythonModel::evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const {
      if (this->user_model.constraint_functions.has_value()) {
         const py::object user_data = this->user_model.user_data.has_value() ? *this->user_model.user_data : py::cast(nullptr);
         const int32_t return_code = (*this->user_model.constraint_functions)(static_cast<int32_t>(this->number_variables),
            static_cast<int32_t>(this->number_constraints), x, constraints.data(), user_data);
         if (0 < return_code) {
            throw FunctionEvaluationError();
         }
      }
   }

   void PythonModel::evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const {
      if (this->user_model.objective_gradient.has_value()) {
         const py::object user_data = this->user_model.user_data.has_value() ? *this->user_model.user_data : py::cast(nullptr);
         const int32_t return_code = (*this->user_model.objective_gradient)(static_cast<int32_t>(this->number_variables),
            x, gradient.data(), user_data);
         if (0 < return_code) {
            throw GradientEvaluationError();
         }
         for (size_t variable_index: Range(this->number_variables)) {
            gradient[variable_index] *= this->optimization_sense;
         }
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
      if (this->user_model.constraint_jacobian.has_value()) {
         const py::object user_data = this->user_model.user_data.has_value() ? *this->user_model.user_data : py::cast(nullptr);
         const int32_t return_code = (*this->user_model.constraint_jacobian)(static_cast<int32_t>(this->number_variables),
            static_cast<int32_t>(this->number_jacobian_nonzeros()), x, jacobian_values, user_data);
         if (0 < return_code) {
            throw GradientEvaluationError();
         }
      }
   }

   void PythonModel::evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         double* hessian_values) const {
      if (this->user_model.lagrangian_hessian.has_value()) {
         objective_multiplier *= this->optimization_sense;
         // if the model has a different sign convention for the Lagrangian than Uno, flip the signs of the multipliers
         if (this->user_model.lagrangian_sign_convention == UNO_MULTIPLIER_POSITIVE) {
            const_cast<Vector<double>&>(multipliers).scale(-1.);
         }
         const py::object user_data = this->user_model.user_data.has_value() ? *this->user_model.user_data : py::cast(nullptr);
         const int32_t return_code = (*this->user_model.lagrangian_hessian)(static_cast<int32_t>(this->number_variables),
            static_cast<int32_t>(this->number_constraints), static_cast<int32_t>(this->number_hessian_nonzeros()), x,
            objective_multiplier, multipliers, hessian_values, user_data);
         // flip the signs of the multipliers back
         if (this->user_model.lagrangian_sign_convention == UNO_MULTIPLIER_POSITIVE) {
            const_cast<Vector<double>&>(multipliers).scale(-1.);
         }
         if (0 < return_code) {
            throw HessianEvaluationError();
         }
      }
      else {
         throw std::runtime_error("evaluate_lagrangian_hessian not implemented");
      }
   }

   void PythonModel::compute_jacobian_vector_product(const double* x, const double* vector, double* result) const {
      if (this->user_model.jacobian_operator.has_value()) {
         const py::object user_data = this->user_model.user_data.has_value() ? *this->user_model.user_data : py::cast(nullptr);
         const int32_t return_code = (*this->user_model.jacobian_operator)(static_cast<int32_t>(this->number_variables),
            static_cast<int32_t>(this->number_constraints), x, true, vector, result, user_data);
         if (0 < return_code) {
            throw GradientEvaluationError();
         }
      }
      else {
         throw std::runtime_error("compute_jacobian_vector_product not implemented");
      }
   }

   void PythonModel::compute_jacobian_transposed_vector_product(const double* x, const double* vector, double* result) const {
      if (this->user_model.jacobian_transposed_operator.has_value()) {
         const py::object user_data = this->user_model.user_data.has_value() ? *this->user_model.user_data : py::cast(nullptr);
         const int32_t return_code = (*this->user_model.jacobian_transposed_operator)(static_cast<int32_t>(this->number_variables),
            static_cast<int32_t>(this->number_constraints), x, true, vector, result, user_data);
         if (0 < return_code) {
            throw GradientEvaluationError();
         }
      }
      else {
         throw std::runtime_error("compute_jacobian_transposed_vector_product not implemented");
      }
   }

   void PythonModel::compute_hessian_vector_product(const double* x, const double* vector, double objective_multiplier,
         const Vector<double>& multipliers, double* result) const {
      if (this->user_model.lagrangian_hessian_operator.has_value()) {
         objective_multiplier *= this->optimization_sense;
         // if the model has a different sign convention for the Lagrangian than Uno, flip the signs of the multipliers
         if (this->user_model.lagrangian_sign_convention == UNO_MULTIPLIER_POSITIVE) {
            const_cast<Vector<double>&>(multipliers).scale(-1.);
         }
         const py::object user_data = this->user_model.user_data.has_value() ? *this->user_model.user_data : py::cast(nullptr);
         const int32_t return_code = (*this->user_model.lagrangian_hessian_operator)(static_cast<int32_t>(this->number_variables),
            static_cast<int32_t>(this->number_constraints), x, true, objective_multiplier, multipliers, vector, result, user_data);
         // flip the signs of the multipliers back
         if (this->user_model.lagrangian_sign_convention == UNO_MULTIPLIER_POSITIVE) {
            const_cast<Vector<double>&>(multipliers).scale(-1.);
         }
         if (0 < return_code) {
            throw HessianEvaluationError();
         }
      }
      else {
         throw std::runtime_error("compute_hessian_vector_product not implemented");
      }
   }

   double PythonModel::variable_lower_bound(size_t variable_index) const {
      if (this->user_model.variables_lower_bounds.has_value()) {
         return (*this->user_model.variables_lower_bounds)[variable_index];
      }
      return -INF<double>;
   }

   double PythonModel::variable_upper_bound(size_t variable_index) const {
      if (this->user_model.variables_upper_bounds.has_value()) {
         return (*this->user_model.variables_upper_bounds)[variable_index];
      }
      return INF<double>;
   }

   const SparseVector<size_t>& PythonModel::get_slacks() const {
      return this->slacks;
   }

   const Vector<size_t>& PythonModel::get_fixed_variables() const {
      return this->fixed_variables;
   }

   double PythonModel::constraint_lower_bound(size_t constraint_index) const {
      if (this->user_model.constraints_lower_bounds.has_value()) {
         return (*this->user_model.constraints_lower_bounds)[constraint_index];
      }
      return -INF<double>;
   }

   double PythonModel::constraint_upper_bound(size_t constraint_index) const {
      if (this->user_model.constraints_upper_bounds.has_value()) {
         return (*this->user_model.constraints_upper_bounds)[constraint_index];
      }
      return INF<double>;
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
         if (this->user_model.lagrangian_sign_convention == UNO_MULTIPLIER_POSITIVE) {
            multipliers.scale(-1.);
         }
      }
      else {
         multipliers.fill(0.);
      }
   }

   void PythonModel::postprocess_solution(Iterate& iterate) const {
      // flip the signs of the multipliers, depending on what the sign convention of the Lagrangian is, and whether
      // we maximize
      iterate.multipliers.constraints *= -this->user_model.lagrangian_sign_convention * this->optimization_sense;
      iterate.multipliers.lower_bounds *= -this->user_model.lagrangian_sign_convention * this->optimization_sense;
      iterate.multipliers.upper_bounds *= -this->user_model.lagrangian_sign_convention * this->optimization_sense;
      iterate.evaluations.objective *= this->optimization_sense;
   }

   size_t PythonModel::number_jacobian_nonzeros() const {
      return static_cast<size_t>(this->user_model.number_jacobian_nonzeros);
   }

   size_t PythonModel::number_hessian_nonzeros() const {
      return static_cast<size_t>(this->user_model.number_hessian_nonzeros);
   }
} // namespace