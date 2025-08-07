// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "HomogeneousEqualityConstrainedModel.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/Concatenation.hpp"
#include "symbolic/CollectionAdapter.hpp"
#include "symbolic/Range.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   // Transform the problem into an equality-constrained problem with constraints c(x) = 0. This implies:
   // - inequality constraints get a slack
   // - equality constraints are shifted by their RHS
   HomogeneousEqualityConstrainedModel::HomogeneousEqualityConstrainedModel(std::unique_ptr<Model> original_model):
         Model(original_model->name + " -> equality constrained", original_model->number_variables +
            original_model->get_inequality_constraints().size(), original_model->number_constraints, original_model->objective_sign),
         // transfer ownership of the pointer
         model(std::move(original_model)),
         constraint_index_of_inequality_index(this->model->get_inequality_constraints().size()),
         slack_index_of_constraint_index(this->model->number_constraints),
         // all constraints are equality constraints
         equality_constraints(Range(this->number_constraints)),
         inequality_constraints(Range(0)),
         slacks(this->model->get_inequality_constraints().size()),
         lower_bounded_slacks(this->slacks.size()),
         upper_bounded_slacks(this->slacks.size()),
         lower_bounded_variables(concatenate(this->model->get_lower_bounded_variables(), adapt(this->lower_bounded_slacks))),
         upper_bounded_variables(concatenate(this->model->get_upper_bounded_variables(), adapt(this->upper_bounded_slacks))),
         single_lower_bounded_variables(concatenate(this->model->get_single_lower_bounded_variables(), adapt(this->single_lower_bounded_slacks))),
         single_upper_bounded_variables(concatenate(this->model->get_single_upper_bounded_variables(), adapt(this->single_upper_bounded_slacks))){
      // register the inequality constraint of each slack
      size_t inequality_index = 0;
      for (const size_t constraint_index: this->model->get_inequality_constraints()) {
         const size_t slack_variable_index = this->model->number_variables + inequality_index;
         this->constraint_index_of_inequality_index[inequality_index] = constraint_index;
         this->slack_index_of_constraint_index[constraint_index] = slack_variable_index;
         this->slacks.insert(constraint_index, slack_variable_index);
         if (is_finite(this->model->constraint_lower_bound(constraint_index))) {
            this->lower_bounded_slacks.emplace_back(slack_variable_index);
            if (!is_finite(this->model->constraint_upper_bound(constraint_index))) {
               this->single_lower_bounded_slacks.emplace_back(slack_variable_index);
            }
         }
         if (is_finite(this->model->constraint_upper_bound(constraint_index))) {
            this->upper_bounded_slacks.emplace_back(slack_variable_index);
            if (!is_finite(this->model->constraint_lower_bound(constraint_index))) {
               this->single_upper_bounded_slacks.emplace_back(slack_variable_index);
            }
         }
         inequality_index++;
      }
   }

   double HomogeneousEqualityConstrainedModel::evaluate_objective(const Vector<double>& x) const {
      return this->model->evaluate_objective(x);
   }

   void HomogeneousEqualityConstrainedModel::evaluate_constraints(const Vector<double>& x, std::vector<double>& constraints) const {
      this->model->evaluate_constraints(x, constraints);
      // inequality constraints: add the slacks
      for (const auto [constraint_index, slack_index]: this->get_slacks()) {
         constraints[constraint_index] -= x[slack_index];
      }

      // equality constraints: make sure they are homogeneous (c(x) = 0)
      for (const size_t constraint_index: this->model->get_equality_constraints()) {
         constraints[constraint_index] -= this->model->constraint_lower_bound(constraint_index);
      }
   }

   void HomogeneousEqualityConstrainedModel::evaluate_objective_gradient(const Vector<double>& x, Vector<double>& gradient) const {
      this->model->evaluate_objective_gradient(x, gradient);
   }

   void HomogeneousEqualityConstrainedModel::compute_jacobian_sparsity(size_t* row_indices, size_t* column_indices,
         size_t solver_indexing) const {
      this->model->compute_jacobian_sparsity(row_indices, column_indices, solver_indexing);

      // add the slack contributions
      size_t nonzero_index = this->model->number_jacobian_nonzeros();
      for (const auto [constraint_index, slack_index]: this->get_slacks()) {
         row_indices[nonzero_index] = constraint_index + solver_indexing;
         column_indices[nonzero_index] = slack_index + solver_indexing;
         ++nonzero_index;
      }
   }

   void HomogeneousEqualityConstrainedModel::compute_hessian_sparsity(size_t* row_indices, size_t* column_indices, size_t solver_indexing) const {
      this->model->compute_hessian_sparsity(row_indices, column_indices, solver_indexing);
   }

   void HomogeneousEqualityConstrainedModel::evaluate_constraint_jacobian(const Vector<double>& x, double* jacobian_values) const {
      this->model->evaluate_constraint_jacobian(x, jacobian_values);

      // add the slack contributions
      size_t nonzero_index = this->model->number_jacobian_nonzeros();
      for (const auto _: this->get_slacks()) {
         jacobian_values[nonzero_index] = -1.;
         ++nonzero_index;
      }
   }

   void HomogeneousEqualityConstrainedModel::evaluate_lagrangian_hessian(const Vector<double>& x, double objective_multiplier,
         const Vector<double>& multipliers, Vector<double>& hessian_values) const {
      this->model->evaluate_lagrangian_hessian(x, objective_multiplier, multipliers, hessian_values);
   }

   void HomogeneousEqualityConstrainedModel::compute_hessian_vector_product(const double* vector, double objective_multiplier,
         const Vector<double>& multipliers, double* result) const {
      this->model->compute_hessian_vector_product(vector, objective_multiplier, multipliers, result);
   }

   double HomogeneousEqualityConstrainedModel::variable_lower_bound(size_t variable_index) const {
      if (variable_index < this->model->number_variables) { // original variable
         return this->model->variable_lower_bound(variable_index);
      }
      else { // slack variable
         const size_t slack_index = variable_index - this->model->number_variables;
         const size_t constraint_index = this->constraint_index_of_inequality_index[slack_index];
         return this->model->constraint_lower_bound(constraint_index);
      }
   }

   double HomogeneousEqualityConstrainedModel::variable_upper_bound(size_t variable_index) const {
      if (variable_index < this->model->number_variables) { // original variable
         return this->model->variable_upper_bound(variable_index);
      }
      else { // slack variable
         const size_t inequality_index = variable_index - this->model->number_variables;
         const size_t constraint_index = this->constraint_index_of_inequality_index[inequality_index];
         return this->model->constraint_upper_bound(constraint_index);
      }
   }

   const Collection<size_t>& HomogeneousEqualityConstrainedModel::get_lower_bounded_variables() const {
      return this->lower_bounded_variables;
   }

   const Collection<size_t>& HomogeneousEqualityConstrainedModel::get_upper_bounded_variables() const {
      return this->upper_bounded_variables;
   }

   const SparseVector<size_t>& HomogeneousEqualityConstrainedModel::get_slacks() const {
      return this->slacks;
   }

   const Collection<size_t>& HomogeneousEqualityConstrainedModel::get_single_lower_bounded_variables() const {
      return this->single_lower_bounded_variables;
   }

   const Collection<size_t>& HomogeneousEqualityConstrainedModel::get_single_upper_bounded_variables() const {
      return this->single_upper_bounded_variables;
   }

   double HomogeneousEqualityConstrainedModel::constraint_lower_bound(size_t /*constraint_index*/) const {
      return 0.; // c(x) = 0
   }

   double HomogeneousEqualityConstrainedModel::constraint_upper_bound(size_t /*constraint_index*/) const {
      return 0.;
   }

   const Collection<size_t>& HomogeneousEqualityConstrainedModel::get_equality_constraints() const {
      return this->equality_constraints;
   }

   const Collection<size_t>& HomogeneousEqualityConstrainedModel::get_inequality_constraints() const {
      return this->inequality_constraints;
   }

   const Collection<size_t>& HomogeneousEqualityConstrainedModel::get_linear_constraints() const {
      return this->model->get_linear_constraints();
   }

   const Vector<size_t>& HomogeneousEqualityConstrainedModel::get_fixed_variables() const {
      return this->model->get_fixed_variables();
   }

   void HomogeneousEqualityConstrainedModel::initial_primal_point(Vector<double>& x) const {
      this->model->initial_primal_point(x);
      // set the slacks
      for (const auto [_, slack_index]: this->get_slacks()) {
         x[slack_index] = 0.;
      }
   }

   void HomogeneousEqualityConstrainedModel::initial_dual_point(Vector<double>& multipliers) const {
      this->model->initial_dual_point(multipliers);
   }

   void HomogeneousEqualityConstrainedModel::postprocess_solution(Iterate& iterate, IterateStatus termination_status) const {
      // discard the slacks
      iterate.number_variables = this->model->number_variables;
      this->model->postprocess_solution(iterate, termination_status);
   }

   size_t HomogeneousEqualityConstrainedModel::number_jacobian_nonzeros() const {
      return this->model->number_jacobian_nonzeros() + this->slacks.size();
   }

   size_t HomogeneousEqualityConstrainedModel::number_hessian_nonzeros() const {
      return this->model->number_hessian_nonzeros();
   }
} // namespace