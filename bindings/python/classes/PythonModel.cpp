// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "PythonModel.hpp"
#include "symbolic/Concatenation.hpp"
#include "Uno.hpp"

namespace uno {
   PythonModel::PythonModel(const std::string& file_name, size_t number_variables, size_t number_constraints,
            double objective_sign, const objective_function_type& objective, const constraint_functions_type& constraints,
            const objective_gradient_type& evaluate_objective_gradient, const jacobian_type& evaluate_jacobian,
            const lagrangian_hessian_type& evaluate_lagrangian_hessian, size_t number_jacobian_nonzeros,
            size_t number_hessian_nonzeros, const std::vector<double>& variables_lower_bounds,
            const std::vector<double>& variables_upper_bounds, const std::vector<double>& constraints_lower_bounds,
            const std::vector<double>& constraints_upper_bounds, const std::vector<double>& primal_initial_point,
            const std::vector<double>& dual_initial_point) :
         Model(file_name, number_variables, number_constraints, objective_sign),
         // functions
         objective(objective), constraints(constraints), objective_gradient(evaluate_objective_gradient),
         jacobian(evaluate_jacobian), hessian(evaluate_lagrangian_hessian),
         // sparsity
         jacobian_nnz(number_jacobian_nonzeros),
         hessian_nnz(number_hessian_nonzeros),
         // bounds
         variables_lower_bounds(variables_lower_bounds), variables_upper_bounds(variables_upper_bounds),
         constraints_lower_bounds(constraints_lower_bounds), constraints_upper_bounds(constraints_upper_bounds),
         // initial point
         primal_initial_point(primal_initial_point), dual_initial_point(dual_initial_point),
         // additional information on the variables and the constraints
         constraint_type(this->number_constraints),
         equality_constraints_collection(this->equality_constraints),
         inequality_constraints_collection(this->inequality_constraints) {
      // variables
      this->fixed_variables.reserve(this->number_variables);
      //Model::partition_variables(this->fixed_variables);

      // constraints
      this->equality_constraints.reserve(this->number_constraints);
      this->inequality_constraints.reserve(this->number_constraints);
      //Model::partition_constraints(this->equality_constraints, this->inequality_constraints);
   }

   bool PythonModel::has_implicit_hessian_representation() const {
      return false;
   }

   bool PythonModel::has_explicit_hessian_representation() const {
      return true;
   }

   double PythonModel::evaluate_objective(const Vector<double>& x) const {
      return this->objective(wrap_pointer(const_cast<Vector<double>*>(&x)));
   }

   // sparse gradient
   void PythonModel::evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const {
      this->objective_gradient(wrap_pointer(const_cast<Vector<double>*>(&x)), wrap_pointer(&gradient));
   }

   void PythonModel::evaluate_constraints(const Vector<double>& x, std::vector<double>& constraints) const {
      this->constraints(wrap_pointer(const_cast<Vector<double>*>(&x)), wrap_pointer(&constraints));
   }

   void PythonModel::evaluate_constraint_jacobian(const Vector<double>& x, double* jacobian_values) const {
      this->jacobian(wrap_pointer(const_cast<Vector<double>*>(&x)), wrap_pointer(jacobian_values));
   }

   void PythonModel::evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         double* hessian_values) const {
      this->hessian(wrap_pointer(const_cast<Vector<double>*>(&x)), objective_multiplier, wrap_pointer(const_cast<Vector<double>*>(&multipliers)),
         wrap_pointer(hessian_values));
   }

   void PythonModel::compute_hessian_vector_product(const double* /*vector*/, double /*objective_multiplier*/,
         const Vector<double>& /*multipliers*/, double* /*result*/) const {
      throw std::runtime_error("PythonModel::compute_hessian_vector_product not implemented yet");
   }

   double PythonModel::variable_lower_bound(size_t variable_index) const {
      return this->variables_lower_bounds[variable_index];
   }

   double PythonModel::variable_upper_bound(size_t variable_index) const {
      return this->variables_upper_bounds[variable_index];
   }

   const SparseVector<size_t>& PythonModel::get_slacks() const {
      return this->slacks;
   }

   const Vector<size_t>& PythonModel::get_fixed_variables() const {
      return this->fixed_variables;
   }

   double PythonModel::constraint_lower_bound(size_t constraint_index) const {
      return this->constraints_lower_bounds[constraint_index];
   }

   double PythonModel::constraint_upper_bound(size_t constraint_index) const {
      return this->constraints_upper_bounds[constraint_index];
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

   // initial primal point
   void PythonModel::initial_primal_point(Vector<double>& x) const {
      x = this->primal_initial_point;
   }

   // initial dual point
   void PythonModel::initial_dual_point(Vector<double>& multipliers) const {
      multipliers = this->dual_initial_point;
   }

   void PythonModel::postprocess_solution(Iterate& /*iterate*/, IterateStatus /*iterate_status*/) const {
      // do nothing
   }

   size_t PythonModel::number_jacobian_nonzeros() const {
      return this->jacobian_nnz;
   }

   size_t PythonModel::number_hessian_nonzeros() const {
      return this->hessian_nnz;
   }
} // namespace