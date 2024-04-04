// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_HOMOGENEOUSEQUALITYCONSTRAINEDMODEL_H
#define UNO_HOMOGENEOUSEQUALITYCONSTRAINEDMODEL_H

#include "Model.hpp"
#include "tools/Infinity.hpp"
#include "tools/Range.hpp"

// generate an equality-constrained model by:
// - introducing slacks in inequality constraints
// - subtracting the (possibly nonzero) RHS of equality constraints
// all constraints are of the form "c(x) = 0"
class HomogeneousEqualityConstrainedModel: public Model {
public:
   explicit HomogeneousEqualityConstrainedModel(std::unique_ptr<Model> original_model);

   [[nodiscard]] double variable_lower_bound(size_t variable_index) const override;
   [[nodiscard]] double variable_upper_bound(size_t variable_index) const override;
   [[nodiscard]] double constraint_lower_bound(size_t constraint_index) const override;
   [[nodiscard]] double constraint_upper_bound(size_t constraint_index) const override;

   [[nodiscard]] double evaluate_objective(const std::vector<double>& x) const override;
   void evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const override;
   void evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const override;
   void evaluate_constraint_gradient(const std::vector<double>& x, size_t constraint_index, SparseVector<double>& gradient) const override;
   void evaluate_constraint_jacobian(const std::vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const override;
   void evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
         SymmetricMatrix<double>& hessian) const override;

   [[nodiscard]] BoundType get_variable_bound_type(size_t variable_index) const override;
   [[nodiscard]] FunctionType get_constraint_type(size_t constraint_index) const override;
   [[nodiscard]] BoundType get_constraint_bound_type(size_t constraint_index) const override;

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
inline HomogeneousEqualityConstrainedModel::HomogeneousEqualityConstrainedModel(std::unique_ptr<Model> original_model):
      Model(original_model->name + "_equalityconstrained", original_model->number_variables + original_model->inequality_constraints.size(),
            original_model->number_constraints),
      // transfer ownership of the pointer
      original_model(std::move(original_model)),
      inequality_constraint_of_slack(this->original_model->inequality_constraints.size()),
      slack_of_inequality_constraint(this->original_model->number_constraints) {
   // all constraints are now equality constraints
   this->equality_constraints.reserve(this->number_constraints);
   this->inequality_constraints.reserve(0);
   for (size_t constraint_index: Range(this->number_constraints)) {
      this->equality_constraints.push_back(constraint_index);
   }

   // figure out bounded variables
   for (size_t variable_index: this->original_model->lower_bounded_variables) {
      this->lower_bounded_variables.push_back(variable_index);
   }
   for (size_t variable_index: this->original_model->upper_bounded_variables) {
      this->upper_bounded_variables.push_back(variable_index);
   }
   for (size_t variable_index: this->original_model->single_lower_bounded_variables) {
      this->single_lower_bounded_variables.push_back(variable_index);
   }
   for (size_t variable_index: this->original_model->single_upper_bounded_variables) {
      this->single_upper_bounded_variables.push_back(variable_index);
   }
   // register the inequality constraint of each slack
   for (size_t variable_index: Range(this->original_model->inequality_constraints.size())) {
      const size_t constraint_index = this->original_model->inequality_constraints[variable_index];
      const size_t slack_index = variable_index + this->original_model->number_variables;
      this->inequality_constraint_of_slack[variable_index] = constraint_index;
      this->slack_of_inequality_constraint[constraint_index] = variable_index;
      this->slacks.insert(constraint_index, slack_index);
      if (is_finite(this->original_model->constraint_lower_bound(constraint_index))) {
         this->lower_bounded_variables.push_back(slack_index);
         if (not is_finite(this->original_model->constraint_upper_bound(constraint_index))) {
            this->single_lower_bounded_variables.push_back(slack_index);
         }
      }
      if (is_finite(this->original_model->constraint_upper_bound(constraint_index))) {
         this->upper_bounded_variables.push_back(slack_index);
         if (not is_finite(this->original_model->constraint_lower_bound(constraint_index))) {
            this->single_upper_bounded_variables.push_back(slack_index);
         }
      }
   }
}

inline double HomogeneousEqualityConstrainedModel::variable_lower_bound(size_t variable_index) const {
   if (variable_index < this->original_model->number_variables) { // original variable
      return this->original_model->variable_lower_bound(variable_index);
   }
   else { // slack variable
      const size_t slack_index = variable_index - this->original_model->number_variables;
      const size_t constraint_index = this->inequality_constraint_of_slack[slack_index];
      return this->original_model->constraint_lower_bound(constraint_index);
   }
}

inline double HomogeneousEqualityConstrainedModel::variable_upper_bound(size_t variable_index) const {
   if (variable_index < this->original_model->number_variables) { // original variable
      return this->original_model->variable_upper_bound(variable_index);
   }
   else { // slack variable
      const size_t slack_index = variable_index - this->original_model->number_variables;
      const size_t constraint_index = this->inequality_constraint_of_slack[slack_index];
      return this->original_model->constraint_upper_bound(constraint_index);
   }
}

inline double HomogeneousEqualityConstrainedModel::constraint_lower_bound(size_t /*j*/) const {
   // all constraints are of the form "c(x) = 0"
   return 0.;
}

inline double HomogeneousEqualityConstrainedModel::constraint_upper_bound(size_t /*j*/) const {
   // all constraints are of the form "c(x) = 0"
   return 0.;
}

inline double HomogeneousEqualityConstrainedModel::evaluate_objective(const std::vector<double>& x) const {
   return this->original_model->evaluate_objective(x);
}

inline void HomogeneousEqualityConstrainedModel::evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const {
   this->original_model->evaluate_objective_gradient(x, gradient);
}

inline void HomogeneousEqualityConstrainedModel::evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const {
   this->original_model->evaluate_constraints(x, constraints);
   // inequality constraints: add the slacks
   for (size_t variable_index: Range(this->original_model->inequality_constraints.size())) {
      const size_t constraint_index = this->original_model->inequality_constraints[variable_index];
      const size_t slack_index = this->original_model->number_variables + variable_index;
      constraints[constraint_index] -= x[slack_index];
   }
   // equality constraints: make sure they are homogeneous (c(x) = 0)
   for (size_t constraint_index: this->original_model->equality_constraints) {
      constraints[constraint_index] -= this->original_model->constraint_lower_bound(constraint_index);
   }
}

inline void HomogeneousEqualityConstrainedModel::evaluate_constraint_gradient(const std::vector<double>& x, size_t constraint_index, SparseVector<double>& gradient) const {
   this->original_model->evaluate_constraint_gradient(x, constraint_index, gradient);
   // if the original constraint is an inequality, add the slack contribution
   if (this->original_model->get_constraint_bound_type(constraint_index) != EQUAL_BOUNDS) {
      const size_t slack_index = this->slack_of_inequality_constraint[constraint_index];
      gradient.insert(slack_index, -1.);
   }
}

inline void HomogeneousEqualityConstrainedModel::evaluate_constraint_jacobian(const std::vector<double>& x, RectangularMatrix<double>& constraint_jacobian) const {
   this->original_model->evaluate_constraint_jacobian(x, constraint_jacobian);
   // add the slack contributions
   for (size_t variable_index: Range(this->original_model->inequality_constraints.size())) {
      const size_t constraint_index = this->original_model->inequality_constraints[variable_index];
      const size_t slack_index = this->original_model->number_variables + variable_index;
      constraint_jacobian[constraint_index].insert(slack_index, -1.);
   }
}

inline void HomogeneousEqualityConstrainedModel::evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
      SymmetricMatrix<double>& hessian) const {
   this->original_model->evaluate_lagrangian_hessian(x, objective_multiplier, multipliers, hessian);
   // extend the dimension of the Hessian by finalizing the remaining columns (note: the slacks do not enter the Hessian)
   for (size_t constraint_index: Range(this->original_model->number_variables, this->number_variables)) {
      hessian.finalize_column(constraint_index);
   }
}

inline BoundType HomogeneousEqualityConstrainedModel::get_variable_bound_type(size_t variable_index) const {
   if (variable_index < this->original_model->number_variables) { // original variable
      return this->original_model->get_variable_bound_type(variable_index);
   }
   else { // slack variable
      const size_t slack_index = variable_index - this->original_model->number_variables;
      const size_t constraint_index = this->inequality_constraint_of_slack[slack_index];
      return this->original_model->get_constraint_bound_type(constraint_index);
   }
}

inline FunctionType HomogeneousEqualityConstrainedModel::get_constraint_type(size_t constraint_index) const {
   return this->original_model->get_constraint_type(constraint_index);
}

inline BoundType HomogeneousEqualityConstrainedModel::get_constraint_bound_type(size_t /*j*/) const {
   // all constraints are of the form "c(x) = 0"
   return EQUAL_BOUNDS;
}

inline size_t HomogeneousEqualityConstrainedModel::get_number_objective_gradient_nonzeros() const {
   return this->original_model->get_number_objective_gradient_nonzeros();
}

inline size_t HomogeneousEqualityConstrainedModel::get_number_jacobian_nonzeros() const {
   return this->original_model->get_number_jacobian_nonzeros();
}

inline size_t HomogeneousEqualityConstrainedModel::get_number_hessian_nonzeros() const {
   return this->original_model->get_number_hessian_nonzeros();
}

inline void HomogeneousEqualityConstrainedModel::get_initial_primal_point(std::vector<double>& x) const {
   this->original_model->get_initial_primal_point(x);
   // set the slacks
   for (size_t variable_index: Range(this->original_model->inequality_constraints.size())) {
      const size_t slack_index = this->original_model->number_variables + variable_index;
      x[slack_index] = 0.;
   }
}

inline void HomogeneousEqualityConstrainedModel::get_initial_dual_point(std::vector<double>& multipliers) const {
   this->original_model->get_initial_dual_point(multipliers);
}

inline void HomogeneousEqualityConstrainedModel::postprocess_solution(Iterate& iterate, TerminationStatus termination_status) const {
   this->original_model->postprocess_solution(iterate, termination_status);

   // discard the slacks
   iterate.number_variables = this->original_model->number_variables;
}

inline const std::vector<size_t>& HomogeneousEqualityConstrainedModel::get_linear_constraints() const {
   return this->original_model->get_linear_constraints();
}

#endif // UNO_HOMOGENEOUSEQUALITYCONSTRAINEDMODEL_H
