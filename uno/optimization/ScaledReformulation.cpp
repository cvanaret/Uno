#include "ScaledReformulation.hpp"

ScaledReformulation::ScaledReformulation(const Problem& original_problem, const Scaling& scaling):
      Problem(original_problem.name + "_scaled", // name
            original_problem.number_variables, // number of variables
            original_problem.number_constraints, // number of constraints
            original_problem.problem_type), // problem type
      original_problem(original_problem),
      scaling(scaling) {
   // check the scaling factors
   assert(0 <= this->scaling.get_objective_scaling() && "Objective scaling failed.");
   for (size_t j = 0; j < this->number_constraints; j++) {
      assert(0 <= this->scaling.get_constraint_scaling(j) && "Constraint scaling failed.");
   }
   // the constraint distribution is the same as in the original problem
   this->equality_constraints = this->original_problem.equality_constraints;
   this->inequality_constraints = this->original_problem.inequality_constraints;
   this->linear_constraints = this->original_problem.linear_constraints;
}

inline double ScaledReformulation::get_variable_lower_bound(size_t i) const {
   return this->original_problem.get_variable_lower_bound(i);
}

inline double ScaledReformulation::get_variable_upper_bound(size_t i) const {
   return this->original_problem.get_variable_upper_bound(i);
}

inline double ScaledReformulation::get_constraint_lower_bound(size_t j) const {
   const double lb = this->original_problem.get_constraint_lower_bound(j);
   // scale
   return this->scaling.get_constraint_scaling(j)*lb;
}

inline double ScaledReformulation::get_constraint_upper_bound(size_t j) const {
   const double ub = this->original_problem.get_constraint_upper_bound(j);
   // scale
   return this->scaling.get_constraint_scaling(j)*ub;
}

inline double ScaledReformulation::evaluate_objective(const std::vector<double>& x) const {
   const double objective = this->original_problem.evaluate_objective(x);
   // scale
   return this->scaling.get_objective_scaling()*objective;
}

inline void ScaledReformulation::evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const {
   this->original_problem.evaluate_objective_gradient(x, gradient);
   // scale
   scale(gradient, this->scaling.get_objective_scaling());
}

inline void ScaledReformulation::evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const {
   this->original_problem.evaluate_constraints(x, constraints);
   // scale
   for (size_t j = 0; j < this->number_constraints; j++) {
      constraints[j] *= this->scaling.get_constraint_scaling(j);
   }
}

inline void ScaledReformulation::evaluate_constraint_gradient(const std::vector<double>& x, size_t j, SparseVector<double>& gradient) const {
   this->original_problem.evaluate_constraint_gradient(x, j, gradient);
   // scale
   scale(gradient, this->scaling.get_constraint_scaling(j));
}

inline void ScaledReformulation::evaluate_constraint_jacobian(const std::vector<double>& x, std::vector<SparseVector<double>>& constraint_jacobian) const {
   for (size_t j = 0; j < this->number_constraints; j++) {
      // scaled evaluation
      this->evaluate_constraint_gradient(x, j, constraint_jacobian[j]);
   }
}

inline void ScaledReformulation::evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier,
      const std::vector<double>& multipliers, SymmetricMatrix& hessian) const {
   // scale the objective and constraint multipliers
   const double scaled_objective_multiplier = objective_multiplier*this->scaling.get_objective_scaling();
   // TODO preallocate this vector
   static std::vector<double> scaled_multipliers(this->number_constraints);
   for (size_t j = 0; j < this->number_constraints; j++) {
      scaled_multipliers[j] = scaling.get_constraint_scaling(j)*multipliers[j];
   }
   this->original_problem.evaluate_lagrangian_hessian(x, scaled_objective_multiplier, scaled_multipliers, hessian);
}

inline ConstraintType ScaledReformulation::get_variable_status(size_t i) const {
   return this->original_problem.get_variable_status(i);
}

inline FunctionType ScaledReformulation::get_constraint_type(size_t j) const {
   return this->original_problem.get_constraint_type(j);
}

inline ConstraintType ScaledReformulation::get_constraint_status(size_t j) const {
   return this->original_problem.get_constraint_status(j);
}

inline size_t ScaledReformulation::get_hessian_maximum_number_nonzeros() const {
   return this->original_problem.get_hessian_maximum_number_nonzeros();
}

inline void ScaledReformulation::get_initial_primal_point(std::vector<double>& x) const {
   this->original_problem.get_initial_primal_point(x);
}

inline void ScaledReformulation::get_initial_dual_point(std::vector<double>& multipliers) const {
   this->original_problem.get_initial_dual_point(multipliers);
}