// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "PythonModel.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/Concatenation.hpp"
#include "Uno.hpp"

namespace uno {
   PythonModel::PythonModel(const std::string& file_name, size_t number_variables, size_t number_constraints,
            double objective_sign, const objective_function_type& evaluate_objective, const constraint_functions_type& evaluate_constraints,
            const objective_gradient_type& evaluate_objective_gradient, const jacobian_type& evaluate_jacobian,
            const lagrangian_hessian_type& evaluate_lagrangian_hessian, const std::vector<double>& variables_lower_bounds,
            const std::vector<double>& variables_upper_bounds, const std::vector<double>& constraints_lower_bounds,
            const std::vector<double>& constraints_upper_bounds, const std::vector<double>& primal_initial_point,
            const std::vector<double>& dual_initial_point) :
         Model(file_name, number_variables, number_constraints, objective_sign),
         // functions
         objective(evaluate_objective), constraints(evaluate_constraints), objective_gradient(evaluate_objective_gradient),
         jacobian(evaluate_jacobian), hessian(evaluate_lagrangian_hessian),
         // bounds
         variables_lower_bounds(variables_lower_bounds), variables_upper_bounds(variables_upper_bounds),
         constraints_lower_bounds(constraints_lower_bounds), constraints_upper_bounds(constraints_upper_bounds),
         // initial point
         primal_initial_point(primal_initial_point), dual_initial_point(dual_initial_point),
         // additional information on the variables and the constraints
         variable_status(this->number_variables),
         constraint_type(this->number_constraints),
         constraint_status(this->number_constraints),
         linear_constraints_collection(this->linear_constraints),
         equality_constraints_collection(this->equality_constraints),
         inequality_constraints_collection(this->inequality_constraints),
         lower_bounded_variables_collection(this->lower_bounded_variables),
         upper_bounded_variables_collection(this->upper_bounded_variables),
         single_lower_bounded_variables_collection(this->single_lower_bounded_variables),
         single_upper_bounded_variables_collection(this->single_upper_bounded_variables) {
      // variables
      this->lower_bounded_variables.reserve(this->number_variables);
      this->upper_bounded_variables.reserve(this->number_variables);
      this->single_lower_bounded_variables.reserve(this->number_variables);
      this->single_upper_bounded_variables.reserve(this->number_variables);
      this->fixed_variables.reserve(this->number_variables);
      for (size_t variable_index: Range(this->number_variables)) {
         if (variables_lower_bounds[variable_index] == variables_upper_bounds[variable_index]) {
            WARNING << "Variable x" << variable_index << " has identical bounds\n";
            this->fixed_variables.emplace_back(variable_index);
         }
      }
      this->categorize_bounded_variables();

      // constraints
      this->equality_constraints.reserve(this->number_constraints);
      this->inequality_constraints.reserve(this->number_constraints);
      this->linear_constraints.reserve(this->number_constraints);
      Model::determine_bounds_types(this->constraints_lower_bounds, this->constraints_upper_bounds, this->constraint_status);
      Model::partition_constraints(this->equality_constraints, this->inequality_constraints);
   }

   double PythonModel::evaluate_objective(const Vector<double>& x) const {
      return this->objective(wrap_pointer(const_cast<Vector<double>*>(&x)));
   }

   // sparse gradient
   void PythonModel::evaluate_objective_gradient(const Vector<double>& x, SparseVector<double>& gradient) const {
      this->objective_gradient(wrap_pointer(const_cast<Vector<double>*>(&x)), wrap_pointer(&gradient));
   }

   void PythonModel::evaluate_constraints(const Vector<double>& x, Vector<double>& constraints) const {
      this->constraints(wrap_pointer(const_cast<Vector<double>*>(&x)), wrap_pointer(&constraints));
   }

   // sparse gradient
   void PythonModel::evaluate_constraint_gradient(const Vector<double>& x, size_t constraint_index, SparseVector<double>& gradient) const {
      throw std::runtime_error("PythonModel::evaluate_constraint_gradient not implemented");
   }

   void PythonModel::evaluate_constraint_jacobian(const Vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const {
      this->jacobian(wrap_pointer(const_cast<Vector<double>*>(&x)), wrap_pointer(&constraint_jacobian));
   }

   void PythonModel::evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier, const Vector<double>& multipliers,
         SymmetricMatrix<size_t, double>& hessian) const {
      this->hessian(wrap_pointer(const_cast<Vector<double>*>(&x)), objective_multiplier, wrap_pointer(const_cast<Vector<double>*>(&multipliers)),
         wrap_pointer(&hessian));
   }

   void PythonModel::compute_hessian_vector_product(const double* vector, double objective_multiplier, const Vector<double>& multipliers,
         double* result) const {
      throw std::runtime_error("PythonModel::compute_hessian_vector_product not implemented yet");
   }

   double PythonModel::variable_lower_bound(size_t variable_index) const {
      return this->variables_lower_bounds[variable_index];
   }

   double PythonModel::variable_upper_bound(size_t variable_index) const {
      return this->variables_upper_bounds[variable_index];
   }

   BoundType PythonModel::get_variable_bound_type(size_t variable_index) const {
      const double lb = this->variable_lower_bound(variable_index);
      const double ub = this->variable_upper_bound(variable_index);
      const bool is_lb_finite = is_finite(lb);
      const bool is_ub_finite = is_finite(ub);
      if (is_lb_finite && is_ub_finite) {
         if (lb == ub) {
            return EQUAL_BOUNDS;
         }
         else {
            return BOUNDED_BOTH_SIDES;
         }
      }
      else if (is_lb_finite) {
         return BOUNDED_LOWER;
      }
      else if (is_ub_finite) {
         return BOUNDED_UPPER;
      }
      return UNBOUNDED;
   }

   const Collection<size_t>& PythonModel::get_lower_bounded_variables() const {
      return this->lower_bounded_variables_collection;
   }

   const Collection<size_t>& PythonModel::get_upper_bounded_variables() const {
      return this->upper_bounded_variables_collection;
   }

   const SparseVector<size_t>& PythonModel::get_slacks() const {
      return this->slacks;
   }

   const Collection<size_t>& PythonModel::get_single_lower_bounded_variables() const {
      return this->single_lower_bounded_variables_collection;
   }

   const Collection<size_t>& PythonModel::get_single_upper_bounded_variables() const {
      return this->single_upper_bounded_variables_collection;
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

   FunctionType PythonModel::get_constraint_type(size_t constraint_index) const {
      return this->constraint_type[constraint_index];
   }

   BoundType PythonModel::get_constraint_bound_type(size_t constraint_index) const {
      return this->constraint_status[constraint_index];
   }

   const Collection<size_t>& PythonModel::get_equality_constraints() const {
      return this->equality_constraints_collection;
   }

   const Collection<size_t>& PythonModel::get_inequality_constraints() const {
      return this->inequality_constraints_collection;
   }

   const Collection<size_t>& PythonModel::get_linear_constraints() const {
      return this->linear_constraints_collection;
   }

   // initial primal point
   void PythonModel::initial_primal_point(Vector<double>& x) const {
      x = this->primal_initial_point;
   }

   // initial dual point
   void PythonModel::initial_dual_point(Vector<double>& multipliers) const {
      multipliers = this->dual_initial_point;
   }

   void PythonModel::postprocess_solution(Iterate& iterate, IterateStatus iterate_status) const {
      // do nothing
   }

   size_t PythonModel::number_objective_gradient_nonzeros() const {
      return 2;
   }

   size_t PythonModel::number_jacobian_nonzeros() const {
      return 4;
   }

   size_t PythonModel::number_hessian_nonzeros() const {
      return 3;
   }

   void PythonModel::categorize_bounded_variables() {
      Model::determine_bounds_types(this->variables_lower_bounds, this->variables_upper_bounds, this->variable_status);
      // figure out the bounded variables
      for (size_t variable_index: Range(this->number_variables)) {
         const BoundType status = this->get_variable_bound_type(variable_index);
         if (status == BOUNDED_LOWER || status == BOUNDED_BOTH_SIDES) {
            this->lower_bounded_variables.emplace_back(variable_index);
            if (status == BOUNDED_LOWER) {
               this->single_lower_bounded_variables.emplace_back(variable_index);
            }
         }
         if (status == BOUNDED_UPPER || status == BOUNDED_BOTH_SIDES) {
            this->upper_bounded_variables.emplace_back(variable_index);
            if (status == BOUNDED_UPPER) {
               this->single_upper_bounded_variables.emplace_back(variable_index);
            }
         }
      }
   }
} // namespace