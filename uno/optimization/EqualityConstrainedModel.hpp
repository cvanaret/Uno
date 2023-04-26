// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_EQUALITYCONSTRAINEDMODEL_H
#define UNO_EQUALITYCONSTRAINEDMODEL_H

#include "Model.hpp"
#include "tools/Infinity.hpp"
#include "tools/Range.hpp"

// all constraints are of the form "c(x) = 0"
class EqualityConstrainedModel: public Model {
public:
   explicit EqualityConstrainedModel(std::unique_ptr<Model> original_model);

   [[nodiscard]] double get_variable_lower_bound(size_t i) const override;
   [[nodiscard]] double get_variable_upper_bound(size_t i) const override;
   [[nodiscard]] double get_constraint_lower_bound(size_t j) const override;
   [[nodiscard]] double get_constraint_upper_bound(size_t j) const override;

   [[nodiscard]] double evaluate_objective(const std::vector<double>& x) const override;
   void evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const override;
   void evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const override;
   void evaluate_constraint_gradient(const std::vector<double>& x, size_t j, SparseVector<double>& gradient) const override;
   void evaluate_constraint_jacobian(const std::vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const override;
   void evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
         SymmetricMatrix<double>& hessian) const override;

   [[nodiscard]] BoundType get_variable_bound_type(size_t i) const override;
   [[nodiscard]] FunctionType get_constraint_type(size_t j) const override;
   [[nodiscard]] BoundType get_constraint_bound_type(size_t j) const override;

   [[nodiscard]] size_t get_number_objective_gradient_nonzeros() const override;
   [[nodiscard]] size_t get_number_jacobian_nonzeros() const override;
   [[nodiscard]] size_t get_number_hessian_nonzeros() const override;

   void get_initial_primal_point(std::vector<double>& x) const override;
   void get_initial_dual_point(std::vector<double>& multipliers) const override;
   void postprocess_solution(Iterate& iterate, TerminationStatus termination_status) const override;

   [[nodiscard]] const std::vector<size_t>& get_linear_constraints() const override;

protected:
   std::unique_ptr<Model> original_model;
   std::vector<size_t> inequality_constraint_of_slack;
   std::vector<size_t> slack_of_inequality_constraint;
};

// Transform the problem into an equality-constrained problem with constraints c(x) = 0. This implies:
// - inequality constraints get a slack
// - equality constraints are shifted by their RHS
inline EqualityConstrainedModel::EqualityConstrainedModel(std::unique_ptr<Model> original_model):
      Model(original_model->name + "_equalityconstrained", original_model->number_variables + original_model->inequality_constraints.size(),
            original_model->number_constraints, original_model->problem_type),
      // transfer ownership of the pointer
      original_model(std::move(original_model)),
      inequality_constraint_of_slack(this->original_model->inequality_constraints.size()),
      slack_of_inequality_constraint(this->original_model->number_constraints) {
   // all constraints are now equality constraints
   this->equality_constraints.reserve(this->number_constraints);
   this->inequality_constraints.reserve(0);
   for (size_t j: Range(this->number_constraints)) {
      this->equality_constraints.push_back(j);
   }

   // figure out bounded variables
   for (size_t i: this->original_model->lower_bounded_variables) {
      this->lower_bounded_variables.push_back(i);
   }
   for (size_t i: this->original_model->upper_bounded_variables) {
      this->upper_bounded_variables.push_back(i);
   }
   for (size_t i: this->original_model->single_lower_bounded_variables) {
      this->single_lower_bounded_variables.push_back(i);
   }
   for (size_t i: this->original_model->single_upper_bounded_variables) {
      this->single_upper_bounded_variables.push_back(i);
   }
   // register the inequality constraint of each slack
   for (size_t i: Range(this->original_model->inequality_constraints.size())) {
      const size_t j = this->original_model->inequality_constraints[i];
      const size_t slack_index = i + this->original_model->number_variables;
      this->inequality_constraint_of_slack[i] = j;
      this->slack_of_inequality_constraint[j] = i;
      this->slacks.insert(j, slack_index);
      if (is_finite(this->original_model->get_constraint_lower_bound(j))) {
         this->lower_bounded_variables.push_back(slack_index);
         if (not is_finite(this->original_model->get_constraint_upper_bound(j))) {
            this->single_lower_bounded_variables.push_back(slack_index);
         }
      }
      if (is_finite(this->original_model->get_constraint_upper_bound(j))) {
         this->upper_bounded_variables.push_back(slack_index);
         if (not is_finite(this->original_model->get_constraint_lower_bound(j))) {
            this->single_upper_bounded_variables.push_back(slack_index);
         }
      }
   }
}

inline double EqualityConstrainedModel::get_variable_lower_bound(size_t i) const {
   if (i < this->original_model->number_variables) { // original variable
      return this->original_model->get_variable_lower_bound(i);
   }
   else { // slack variable
      const size_t slack_index = i - this->original_model->number_variables;
      const size_t j = this->inequality_constraint_of_slack[slack_index];
      return this->original_model->get_constraint_lower_bound(j);
   }
}

inline double EqualityConstrainedModel::get_variable_upper_bound(size_t i) const {
   if (i < this->original_model->number_variables) { // original variable
      return this->original_model->get_variable_upper_bound(i);
   }
   else { // slack variable
      const size_t slack_index = i - this->original_model->number_variables;
      const size_t j = this->inequality_constraint_of_slack[slack_index];
      return this->original_model->get_constraint_upper_bound(j);
   }
}

inline double EqualityConstrainedModel::get_constraint_lower_bound(size_t /*j*/) const {
   // all constraints are of the form "c(x) = 0"
   return 0.;
}

inline double EqualityConstrainedModel::get_constraint_upper_bound(size_t /*j*/) const {
   // all constraints are of the form "c(x) = 0"
   return 0.;
}

inline double EqualityConstrainedModel::evaluate_objective(const std::vector<double>& x) const {
   return this->original_model->evaluate_objective(x);
}

inline void EqualityConstrainedModel::evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const {
   this->original_model->evaluate_objective_gradient(x, gradient);
}

inline void EqualityConstrainedModel::evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const {
   this->original_model->evaluate_constraints(x, constraints);
   // inequality constraints: add the slacks
   for (size_t i: Range(this->original_model->inequality_constraints.size())) {
      const size_t j = this->original_model->inequality_constraints[i];
      const size_t slack_index = this->original_model->number_variables + i;
      constraints[j] -= x[slack_index];
   }
   // equality constraints: make sure they are "c(x) = 0"
   for (size_t j: this->original_model->equality_constraints) {
      constraints[j] -= this->original_model->get_constraint_lower_bound(j);
   }
}

inline void EqualityConstrainedModel::evaluate_constraint_gradient(const std::vector<double>& x, size_t j, SparseVector<double>& gradient) const {
   this->original_model->evaluate_constraint_gradient(x, j, gradient);
   // if the original constraint is an inequality, add the slack contribution
   if (this->original_model->get_constraint_bound_type(j) != EQUAL_BOUNDS) {
      const size_t slack_index = this->slack_of_inequality_constraint[j];
      gradient.insert(slack_index, -1.);
   }
}

inline void EqualityConstrainedModel::evaluate_constraint_jacobian(const std::vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const {
   this->original_model->evaluate_constraint_jacobian(x, constraint_jacobian);
   // add the slack contributions
   for (size_t i: Range(this->original_model->inequality_constraints.size())) {
      const size_t j = this->original_model->inequality_constraints[i];
      const size_t slack_index = this->original_model->number_variables + i;
      constraint_jacobian[j].insert(slack_index, -1.);
   }
}

inline void EqualityConstrainedModel::evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
      SymmetricMatrix<double>& hessian) const {
   this->original_model->evaluate_lagrangian_hessian(x, objective_multiplier, multipliers, hessian);
   // extend the dimension of the Hessian by finalizing the remaining columns (note: the slacks do not enter the Hessian)
   for (size_t j: Range(this->original_model->number_variables, this->number_variables)) {
      hessian.finalize_column(j);
   }
}

inline BoundType EqualityConstrainedModel::get_variable_bound_type(size_t i) const {
   if (i < this->original_model->number_variables) { // original variable
      return this->original_model->get_variable_bound_type(i);
   }
   else { // slack variable
      const size_t slack_index = i - this->original_model->number_variables;
      const size_t j = this->inequality_constraint_of_slack[slack_index];
      return this->original_model->get_constraint_bound_type(j);
   }
}

inline FunctionType EqualityConstrainedModel::get_constraint_type(size_t j) const {
   return this->original_model->get_constraint_type(j);
}

inline BoundType EqualityConstrainedModel::get_constraint_bound_type(size_t /*j*/) const {
   // all constraints are of the form "c(x) = 0"
   return EQUAL_BOUNDS;
}

inline size_t EqualityConstrainedModel::get_number_objective_gradient_nonzeros() const {
   return this->original_model->get_number_objective_gradient_nonzeros();
}

inline size_t EqualityConstrainedModel::get_number_jacobian_nonzeros() const {
   return this->original_model->get_number_jacobian_nonzeros();
}

inline size_t EqualityConstrainedModel::get_number_hessian_nonzeros() const {
   return this->original_model->get_number_hessian_nonzeros();
}

inline void EqualityConstrainedModel::get_initial_primal_point(std::vector<double>& x) const {
   this->original_model->get_initial_primal_point(x);
   // set the slacks
   for (size_t i: Range(this->original_model->inequality_constraints.size())) {
      const size_t slack_index = this->original_model->number_variables + i;
      x[slack_index] = 0.;
   }
}

inline void EqualityConstrainedModel::get_initial_dual_point(std::vector<double>& multipliers) const {
   this->original_model->get_initial_dual_point(multipliers);
}

inline void EqualityConstrainedModel::postprocess_solution(Iterate& iterate, TerminationStatus termination_status) const {
   this->original_model->postprocess_solution(iterate, termination_status);

   // discard the slacks
   iterate.number_variables = this->original_model->number_variables;
}

inline const std::vector<size_t>& EqualityConstrainedModel::get_linear_constraints() const {
   return this->original_model->get_linear_constraints();
}

#endif // UNO_EQUALITYCONSTRAINEDMODEL_H