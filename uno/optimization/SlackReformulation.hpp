#ifndef UNO_SLACKREFORMULATION_H
#define UNO_SLACKREFORMULATION_H

#include "Problem.hpp"

class SlackReformulation: public Problem {
public:
   explicit SlackReformulation(const Problem& original_problem);

   [[nodiscard]] double get_variable_lower_bound(size_t i) const override;
   [[nodiscard]] double get_variable_upper_bound(size_t i) const override;
   [[nodiscard]] double get_constraint_lower_bound(size_t j) const override;
   [[nodiscard]] double get_constraint_upper_bound(size_t j) const override;

   [[nodiscard]] double evaluate_objective(const std::vector<double>& x) const override;
   void evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const override;
   void evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const override;
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
   const Problem& original_problem;
   std::vector<size_t> inequality_constraint_of_slack;
};

inline SlackReformulation::SlackReformulation(const Problem& original_problem):
      Problem(original_problem.name + "_slacks", // name
            original_problem.number_variables + original_problem.inequality_constraints.size(), // number of variables
            original_problem.number_constraints, // number of constraints
            original_problem.problem_type), // problem type
      original_problem(original_problem),
      inequality_constraint_of_slack(original_problem.inequality_constraints.size()) {
   // all constraints are now equality constraints
   for (size_t j = 0; j < this->number_constraints; j++) {
      this->equality_constraints.insert(j, j);
   }
   // register the inequality constraint of each slack
   this->original_problem.inequality_constraints.for_each([&](size_t j, size_t i) {
      this->inequality_constraint_of_slack[i] = j;
   });
}

inline double SlackReformulation::get_variable_lower_bound(size_t i) const {
   if (i < this->original_problem.number_variables) { // original variable
      return this->original_problem.get_variable_lower_bound(i);
   }
   else { // slack variable
      const size_t slack_index = i - this->original_problem.number_variables;
      const size_t j = this->inequality_constraint_of_slack[slack_index];
      return this->original_problem.get_constraint_lower_bound(j);
   }
}

inline double SlackReformulation::get_variable_upper_bound(size_t i) const {
   if (i < this->original_problem.number_variables) { // original variable
      return this->original_problem.get_variable_upper_bound(i);
   }
   else { // slack variable
      const size_t slack_index = i - this->original_problem.number_variables;
      const size_t j = this->inequality_constraint_of_slack[slack_index];
      return this->original_problem.get_constraint_upper_bound(j);
   }
}

inline double SlackReformulation::get_constraint_lower_bound(size_t /*j*/) const {
   // all constraints are of the form "c(x) = 0"
   return 0.;
}

inline double SlackReformulation::get_constraint_upper_bound(size_t /*j*/) const {
   // all constraints are of the form "c(x) = 0"
   return 0.;
}

inline double SlackReformulation::evaluate_objective(const std::vector<double>& x) const {
   return this->original_problem.evaluate_objective(x);
}

inline void SlackReformulation::evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const {
   this->original_problem.evaluate_objective_gradient(x, gradient);
}

inline void SlackReformulation::evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const {
   this->original_problem.evaluate_constraints(x, constraints);
   // inequality constraints: add the slacks
   this->original_problem.inequality_constraints.for_each([&](size_t j, size_t i) {
      const size_t slack_index = this->original_problem.number_variables + i;
      constraints[j] -= x[slack_index];
   });
   // make sure the equality constraints are "c(x) = 0"
   this->original_problem.equality_constraints.for_each_index([&](size_t j) {
      constraints[j] -= this->original_problem.get_constraint_lower_bound(j);
   });
}

inline void SlackReformulation::evaluate_constraint_jacobian(const std::vector<double>& x, std::vector<SparseVector<double>>& constraint_jacobian) const {
   this->original_problem.evaluate_constraint_jacobian(x, constraint_jacobian);
   // add the slack contributions
   this->original_problem.inequality_constraints.for_each([&](size_t j, size_t i) {
      const size_t slack_index = this->original_problem.number_variables + i;
      constraint_jacobian[j].insert(slack_index, -1.);
   });
}

inline void SlackReformulation::evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
      SymmetricMatrix& hessian) const {
   this->original_problem.evaluate_lagrangian_hessian(x, objective_multiplier, multipliers, hessian);
   // extend the dimension of the Hessian by finalizing the remaining columns (note: the slacks do not enter the Hessian)
   hessian.dimension = this->number_variables;
   for (size_t j = this->original_problem.number_variables; j < this->number_variables; j++) {
      hessian.finalize(j);
   }
}

inline ConstraintType SlackReformulation::get_variable_status(size_t i) const {
   if (i < this->original_problem.number_variables) { // original variable
      return this->original_problem.get_variable_status(i);
   }
   else { // slack variable
      const size_t slack_index = i - this->original_problem.number_variables;
      const size_t j = this->inequality_constraint_of_slack[slack_index];
      return this->original_problem.get_constraint_status(j);
   }
}

inline FunctionType SlackReformulation::get_constraint_type(size_t j) const {
   return this->original_problem.get_constraint_type(j);
}

inline ConstraintType SlackReformulation::get_constraint_status(size_t /*j*/) const {
   // all constraints are of the form "c(x) = 0"
   return EQUAL_BOUNDS;
}

inline size_t SlackReformulation::get_hessian_maximum_number_nonzeros() const {
   return this->original_problem.get_hessian_maximum_number_nonzeros();
}

inline void SlackReformulation::get_initial_primal_point(std::vector<double>& x) const {
   this->original_problem.get_initial_primal_point(x);
   // add the slacks
   this->original_problem.inequality_constraints.for_each_value([&](size_t i) {
      const size_t slack_index = this->original_problem.number_variables + i;
      x[slack_index] = 0.;
   });
}

inline void SlackReformulation::get_initial_dual_point(std::vector<double>& multipliers) const {
   this->original_problem.get_initial_dual_point(multipliers);
}

#endif // UNO_SLACKREFORMULATION_H