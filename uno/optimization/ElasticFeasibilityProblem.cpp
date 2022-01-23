#include <cmath>
#include "ElasticFeasibilityProblem.hpp"
#include "ingredients/constraint_relaxation/ConstraintRelaxationStrategy.hpp"

ElasticFeasibilityProblem::ElasticFeasibilityProblem(const Problem& original_problem, double objective_multiplier,
         double elastic_objective_coefficient, double proximal_coefficient):
      Problem(original_problem.name + "_slacks", // name
            original_problem.number_variables + ConstraintRelaxationStrategy::count_elastic_variables(original_problem), // number of variables
            original_problem.number_constraints, // number of constraints
            original_problem.problem_type), // problem type
      original_problem(original_problem),
      objective_multiplier(objective_multiplier),
      // elastic variables
      elastic_variables(this->number_constraints),
      elastic_objective_coefficient(elastic_objective_coefficient),
      // proximal term
      proximal_coefficient(proximal_coefficient),
      proximal_reference_point(original_problem.number_variables) {
   // register equality and inequality constraints
   this->original_problem.equality_constraints.for_each([&](size_t j, size_t i) {
      this->equality_constraints.insert(j, i);
   });
   this->original_problem.inequality_constraints.for_each([&](size_t j, size_t i) {
      this->inequality_constraints.insert(j, i);
   });
   // generate elastic variables
   ConstraintRelaxationStrategy::generate_elastic_variables(this->original_problem, this->elastic_variables, this->original_problem.number_variables);
}

inline double ElasticFeasibilityProblem::get_variable_lower_bound(size_t i) const {
   if (i < this->original_problem.number_variables) { // original variable
      return this->original_problem.get_variable_lower_bound(i);
   }
   else { // elastic variable in [0, +inf[
      return 0.;
   }
}

inline double ElasticFeasibilityProblem::get_variable_upper_bound(size_t i) const {
   if (i < this->original_problem.number_variables) { // original variable
      return this->original_problem.get_variable_upper_bound(i);
   }
   else { // elastic variable in [0, +inf[
      return std::numeric_limits<double>::infinity();
   }
}

inline double ElasticFeasibilityProblem::get_constraint_lower_bound(size_t j) const {
   return this->original_problem.get_constraint_lower_bound(j);
}

inline double ElasticFeasibilityProblem::get_constraint_upper_bound(size_t j) const {
   return this->original_problem.get_constraint_upper_bound(j);
}

inline double ElasticFeasibilityProblem::compute_elastic_residual(const std::vector<double>& x) const {
   double residual = 0.;
   // l1 residual of the linearized constraints: sum of elastic variables
   auto elastic_contribution = [&](size_t i) {
      residual += x[i];
   };
   this->elastic_variables.positive.for_each_value(elastic_contribution);
   this->elastic_variables.negative.for_each_value(elastic_contribution);
   return residual;
}

inline double ElasticFeasibilityProblem::get_proximal_weight(size_t i) const {
   // weight of each diagonal term of the proximal term
   return std::min(1., 1. / std::abs(this->proximal_reference_point[i]));
}

// return rho*f(x) + coeff*(e^T p + e^T n) + proximal
inline double ElasticFeasibilityProblem::evaluate_objective(const std::vector<double>& x) const {
   // elastic contribution
   double objective = this->elastic_objective_coefficient*this->compute_elastic_residual(x);
   // original objective
   if (this->objective_multiplier != 0.) {
      objective += this->objective_multiplier*this->original_problem.evaluate_objective(x);
   }
   // proximal term
   for (size_t i = 0; i < this->original_problem.number_variables; i++) {
      const double weight = this->get_proximal_weight(i);
      // weighted distance between trial iterate and current iterate
      objective += this->proximal_coefficient * std::pow(weight * (x[i] - this->proximal_reference_point[i]), 2);
   }
   return objective;
}

inline void ElasticFeasibilityProblem::evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const {
   // scale nabla f(x) by rho
   if (this->objective_multiplier != 0.) {
      this->original_problem.evaluate_objective_gradient(x, gradient);
      scale(gradient, this->objective_multiplier);
   }
   else {
      gradient.clear();
   }
   // elastic contribution
   this->elastic_variables.positive.for_each_value([&](size_t elastic_index) {
      gradient.insert(elastic_index, this->elastic_objective_coefficient);
   });
   // proximal term
   for (size_t i = 0; i < this->original_problem.number_variables; i++) {
      const double weight = this->get_proximal_weight(i);
      // measure weighted distance between trial iterate and current iterate
      const double derivative = this->proximal_coefficient * weight * (x[i] - this->proximal_reference_point[i]);
      gradient.insert(i, derivative);
   }
}

inline void ElasticFeasibilityProblem::evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const {
   this->original_problem.evaluate_constraints(x, constraints);
   // add the contribution of the elastics
   this->elastic_variables.positive.for_each([&](size_t j, size_t elastic_index) {
      constraints[j] -= x[elastic_index];
   });
   this->elastic_variables.negative.for_each([&](size_t j, size_t elastic_index) {
      constraints[j] += x[elastic_index];
   });
}

inline void ElasticFeasibilityProblem::evaluate_constraint_jacobian(const std::vector<double>& x, std::vector<SparseVector<double>>& constraint_jacobian) const {
   this->original_problem.evaluate_constraint_jacobian(x, constraint_jacobian);
   // add the contribution of the elastics
   this->elastic_variables.positive.for_each([&](size_t j, size_t elastic_index) {
      constraint_jacobian[j].insert(elastic_index, -1.);
   });
   this->elastic_variables.negative.for_each([&](size_t j, size_t elastic_index) {
      constraint_jacobian[j].insert(elastic_index, 1.);
   });
}

inline void ElasticFeasibilityProblem::evaluate_lagrangian_hessian(const std::vector<double>& x, double /*objective_multiplier*/,
      const std::vector<double>& multipliers, SymmetricMatrix& hessian) const {
   this->original_problem.evaluate_lagrangian_hessian(x, this->objective_multiplier, multipliers, hessian);
   hessian.dimension = this->number_variables;
   // add proximal term for the original variables
   for (size_t i = 0; i < this->original_problem.number_variables; i++) {
      const double distance = std::pow(this->get_proximal_weight(i), 2);
      const double diagonal_term = this->proximal_coefficient*distance;
      hessian.insert(diagonal_term, i, i);
   }
   // extend the dimension of the Hessian by finalizing the remaining columns (note: the elastics do not enter the Hessian)
   for (size_t j = this->original_problem.number_variables; j < this->number_variables; j++) {
      hessian.finalize(j);
   }
}

inline ConstraintType ElasticFeasibilityProblem::get_variable_status(size_t i) const {
   if (i < this->original_problem.number_variables) { // original variable
      return this->original_problem.get_variable_status(i);
   }
   else { // elastic variable in [0, +inf[
      return BOUNDED_LOWER;
   }
}

inline FunctionType ElasticFeasibilityProblem::get_constraint_type(size_t j) const {
   return this->original_problem.get_constraint_type(j);
}

inline ConstraintType ElasticFeasibilityProblem::get_constraint_status(size_t j) const {
   return this->original_problem.get_constraint_status(j);
}

inline size_t ElasticFeasibilityProblem::get_hessian_maximum_number_nonzeros() const {
   return this->original_problem.get_hessian_maximum_number_nonzeros();
}

inline void ElasticFeasibilityProblem::get_initial_primal_point(std::vector<double>& x) const {
   this->original_problem.get_initial_primal_point(x);
   // add the contribution of the elastics
   this->elastic_variables.positive.for_each_value([&](size_t elastic_index) {
      x[elastic_index] = 0.;
   });
   this->elastic_variables.negative.for_each_value([&](size_t elastic_index) {
      x[elastic_index] = 0.;
   });
}

inline void ElasticFeasibilityProblem::get_initial_dual_point(std::vector<double>& multipliers) const {
   this->original_problem.get_initial_dual_point(multipliers);
}

inline void ElasticFeasibilityProblem::set_objective_multiplier(double new_objective_multiplier) {
   // update the objective multiplier
   this->objective_multiplier = new_objective_multiplier;
}

inline void ElasticFeasibilityProblem::set_proximal_coefficient(double new_proximal_coefficient) {
   // update the proximal coefficient
   this->proximal_coefficient = new_proximal_coefficient;
}

inline void ElasticFeasibilityProblem::set_proximal_reference_point(const std::vector<double>& new_proximal_reference_point) {
   // update the proximal reference point
   for (size_t i = 0; i < this->original_problem.number_variables; i++) {
      this->proximal_reference_point[i] = new_proximal_reference_point[i];
   }
}