#ifndef UNO_EQUALITYCONSTRAINEDMODEL_H
#define UNO_EQUALITYCONSTRAINEDMODEL_H

#include "Model.hpp"

class EqualityConstrainedModel: public Model {
public:
   explicit EqualityConstrainedModel(const Model& original_model);

   [[nodiscard]] size_t get_number_original_variables() const override;
   [[nodiscard]] double get_variable_lower_bound(size_t i) const override;
   [[nodiscard]] double get_variable_upper_bound(size_t i) const override;
   [[nodiscard]] double get_constraint_lower_bound(size_t j) const override;
   [[nodiscard]] double get_constraint_upper_bound(size_t j) const override;

   [[nodiscard]] double evaluate_objective(const std::vector<double>& x) const override;
   void evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const override;
   void evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const override;
   void evaluate_constraint_gradient(const std::vector<double>& x, size_t j, SparseVector<double>& gradient) const override;
   void evaluate_constraint_jacobian(const std::vector<double>& x, std::vector<SparseVector<double>>& constraint_jacobian) const override;
   void evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
         SymmetricMatrix& hessian) const override;

   [[nodiscard]] ConstraintType get_variable_status(size_t i) const override;
   [[nodiscard]] FunctionType get_constraint_type(size_t j) const override;
   [[nodiscard]] ConstraintType get_constraint_status(size_t j) const override;
   [[nodiscard]] size_t get_hessian_maximum_number_nonzeros() const override;

   void get_initial_primal_point(std::vector<double>& x) const override;
   void get_initial_dual_point(std::vector<double>& multipliers) const override;

protected:
   const Model& original_model;
   std::vector<size_t> inequality_constraint_of_slack;
};

inline EqualityConstrainedModel::EqualityConstrainedModel(const Model& original_model):
      Model(original_model.name + "_slacks", // name
            original_model.number_variables + original_model.inequality_constraints.size(), // number of variables
            original_model.number_constraints, // number of constraints
            original_model.problem_type), // problem type
      original_model(original_model),
      inequality_constraint_of_slack(original_model.inequality_constraints.size()) {
   // all constraints are now equality constraints
   for (size_t j = 0; j < this->number_constraints; j++) {
      this->equality_constraints.insert(j, j);
   }

   // figure out bounded variables
   for (size_t i: this->original_model.lower_bounded_variables) {
      this->lower_bounded_variables.push_back(i);
   }
   for (size_t i: this->original_model.upper_bounded_variables) {
      this->upper_bounded_variables.push_back(i);
   }
   this->original_model.inequality_constraints.for_each([&](size_t j, size_t i) {
      const size_t slack_index = i + this->original_model.number_variables;
      if (is_finite(this->original_model.get_constraint_lower_bound(j))) {
         this->lower_bounded_variables.push_back(slack_index);
      }
      if (is_finite(this->original_model.get_constraint_upper_bound(j))) {
         this->upper_bounded_variables.push_back(slack_index);
      }
   });

   // register the inequality constraint of each slack
   this->original_model.inequality_constraints.for_each([&](size_t j, size_t i) {
      this->inequality_constraint_of_slack[i] = j;
   });
}

inline size_t EqualityConstrainedModel::get_number_original_variables() const {
   return this->original_model.get_number_original_variables();
}

inline double EqualityConstrainedModel::get_variable_lower_bound(size_t i) const {
   if (i < this->original_model.number_variables) { // original variable
      return this->original_model.get_variable_lower_bound(i);
   }
   else { // slack variable
      const size_t slack_index = i - this->original_model.number_variables;
      const size_t j = this->inequality_constraint_of_slack[slack_index];
      return this->original_model.get_constraint_lower_bound(j);
   }
}

inline double EqualityConstrainedModel::get_variable_upper_bound(size_t i) const {
   if (i < this->original_model.number_variables) { // original variable
      return this->original_model.get_variable_upper_bound(i);
   }
   else { // slack variable
      const size_t slack_index = i - this->original_model.number_variables;
      const size_t j = this->inequality_constraint_of_slack[slack_index];
      return this->original_model.get_constraint_upper_bound(j);
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
   return this->original_model.evaluate_objective(x);
}

inline void EqualityConstrainedModel::evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const {
   this->original_model.evaluate_objective_gradient(x, gradient);
}

inline void EqualityConstrainedModel::evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const {
   this->original_model.evaluate_constraints(x, constraints);
   // inequality constraints: add the slacks
   this->original_model.inequality_constraints.for_each([&](size_t j, size_t i) {
      const size_t slack_index = this->original_model.number_variables + i;
      constraints[j] -= x[slack_index];
   });
   // make sure the equality constraints are "c(x) = 0"
   this->original_model.equality_constraints.for_each_index([&](size_t j) {
      constraints[j] -= this->original_model.get_constraint_lower_bound(j);
   });
}

inline void EqualityConstrainedModel::evaluate_constraint_gradient(const std::vector<double>& x, size_t j, SparseVector<double>& gradient) const {
   this->original_model.evaluate_constraint_gradient(x, j, gradient);
   // add the possible slack contribution
   assert(false && "EqualityConstrainedModel::evaluate_constraint_gradient TODO");
}

inline void EqualityConstrainedModel::evaluate_constraint_jacobian(const std::vector<double>& x, std::vector<SparseVector<double>>& constraint_jacobian) const {
   this->original_model.evaluate_constraint_jacobian(x, constraint_jacobian);
   // add the slack contributions
   this->original_model.inequality_constraints.for_each([&](size_t j, size_t i) {
      const size_t slack_index = this->original_model.number_variables + i;
      constraint_jacobian[j].insert(slack_index, -1.);
   });
}

inline void EqualityConstrainedModel::evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
      SymmetricMatrix& hessian) const {
   this->original_model.evaluate_lagrangian_hessian(x, objective_multiplier, multipliers, hessian);
   // extend the dimension of the Hessian by finalizing the remaining columns (note: the slacks do not enter the Hessian)
   hessian.dimension = this->number_variables;
   for (size_t j = this->original_model.number_variables; j < this->number_variables; j++) {
      hessian.finalize(j);
   }
}

inline ConstraintType EqualityConstrainedModel::get_variable_status(size_t i) const {
   if (i < this->original_model.number_variables) { // original variable
      return this->original_model.get_variable_status(i);
   }
   else { // slack variable
      const size_t slack_index = i - this->original_model.number_variables;
      const size_t j = this->inequality_constraint_of_slack[slack_index];
      return this->original_model.get_constraint_status(j);
   }
}

inline FunctionType EqualityConstrainedModel::get_constraint_type(size_t j) const {
   return this->original_model.get_constraint_type(j);
}

inline ConstraintType EqualityConstrainedModel::get_constraint_status(size_t /*j*/) const {
   // all constraints are of the form "c(x) = 0"
   return EQUAL_BOUNDS;
}

inline size_t EqualityConstrainedModel::get_hessian_maximum_number_nonzeros() const {
   return this->original_model.get_hessian_maximum_number_nonzeros();
}

inline void EqualityConstrainedModel::get_initial_primal_point(std::vector<double>& x) const {
   this->original_model.get_initial_primal_point(x);
   // add the slacks
   this->original_model.inequality_constraints.for_each_value([&](size_t i) {
      const size_t slack_index = this->original_model.number_variables + i;
      x[slack_index] = 0.;
   });
}

inline void EqualityConstrainedModel::get_initial_dual_point(std::vector<double>& multipliers) const {
   this->original_model.get_initial_dual_point(multipliers);
}

#endif // UNO_EQUALITYCONSTRAINEDMODEL_H