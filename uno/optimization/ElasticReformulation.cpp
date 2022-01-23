#include "ElasticReformulation.hpp"
#include "ingredients/constraint_relaxation/ConstraintRelaxationStrategy.hpp"

ElasticReformulation::ElasticReformulation(const Problem& original_problem, double objective_multiplier):
      Problem(original_problem.name + "_slacks", // name
            original_problem.number_variables + original_problem.inequality_constraints.size(), // number of variables
            original_problem.number_constraints, // number of constraints
            original_problem.problem_type), // problem type
      original_problem(original_problem),
      objective_multiplier(objective_multiplier),
      elastic_variables(this->number_variables) {
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

inline double ElasticReformulation::get_variable_lower_bound(size_t i) const {
   if (i < this->original_problem.number_variables) { // original variable
      return this->original_problem.get_variable_lower_bound(i);
   }
   else { // elastic variable in [0, +inf[
      return 0.;
   }
}

inline double ElasticReformulation::get_variable_upper_bound(size_t i) const {
   if (i < this->original_problem.number_variables) { // original variable
      return this->original_problem.get_variable_upper_bound(i);
   }
   else { // elastic variable in [0, +inf[
      return std::numeric_limits<double>::infinity();
   }
}

inline double ElasticReformulation::get_constraint_lower_bound(size_t j) const {
   return this->original_problem.get_constraint_lower_bound(j);
}

inline double ElasticReformulation::get_constraint_upper_bound(size_t j) const {
   return this->original_problem.get_constraint_upper_bound(j);
}

inline double ElasticReformulation::compute_elastic_residual(const std::vector<double>& x) const {
   double residual = 0.;
   // l1 residual of the linearized constraints: sum of elastic variables
   auto elastic_contribution = [&](size_t i) {
      residual += x[i];
   };
   this->elastic_variables.positive.for_each_value(elastic_contribution);
   this->elastic_variables.negative.for_each_value(elastic_contribution);
   return residual;
}

inline double ElasticReformulation::evaluate_objective(const std::vector<double>& x) const {
   // return rho*f(x) + e^T p + e^T n
   double objective = this->compute_elastic_residual(x);
   if (this->objective_multiplier != 0.) {
      objective += this->objective_multiplier*this->original_problem.evaluate_objective(x);
   }
   return objective;
}

inline void ElasticReformulation::evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const {
   // scale nabla f(x) by rho
   if (this->objective_multiplier != 0.) {
      this->original_problem.evaluate_objective_gradient(x, gradient);
      scale(gradient, this->objective_multiplier);
   }
   else {
      gradient.clear();
   }
   // add the contribution of the elastics
   this->elastic_variables.positive.for_each_value([&](size_t elastic_index) {
      gradient.insert(elastic_index, 1.);
   });
}

inline void ElasticReformulation::evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const {
   this->original_problem.evaluate_constraints(x, constraints);
   // add the contribution of the elastics
   this->elastic_variables.positive.for_each([&](size_t j, size_t elastic_index) {
      constraints[j] -= x[elastic_index];
   });
   this->elastic_variables.negative.for_each([&](size_t j, size_t elastic_index) {
      constraints[j] += x[elastic_index];
   });
}

inline void ElasticReformulation::evaluate_constraint_jacobian(const std::vector<double>& x, std::vector<SparseVector<double>>& constraint_jacobian) const {
   this->original_problem.evaluate_constraint_jacobian(x, constraint_jacobian);
   // add the contribution of the elastics
   this->elastic_variables.positive.for_each([&](size_t j, size_t elastic_index) {
      constraint_jacobian[j].insert(elastic_index, -1.);
   });
   this->elastic_variables.negative.for_each([&](size_t j, size_t elastic_index) {
      constraint_jacobian[j].insert(elastic_index, 1.);
   });
}

inline void ElasticReformulation::evaluate_lagrangian_hessian(const std::vector<double>& x, double /*objective_multiplier*/,
      const std::vector<double>& multipliers, SymmetricMatrix& hessian) const {
   this->original_problem.evaluate_lagrangian_hessian(x, this->objective_multiplier, multipliers, hessian);
   // extend the dimension of the Hessian by finalizing the remaining columns (note: the elastics do not enter the Hessian)
   hessian.dimension = this->number_variables;
   for (size_t j = this->original_problem.number_variables; j < this->number_variables; j++) {
      hessian.finalize(j);
   }
}

inline ConstraintType ElasticReformulation::get_variable_status(size_t i) const {
   if (i < this->original_problem.number_variables) { // original variable
      return this->original_problem.get_variable_status(i);
   }
   else { // elastic variable in [0, +inf[
      return BOUNDED_LOWER;
   }
}

inline FunctionType ElasticReformulation::get_constraint_type(size_t j) const {
   return this->original_problem.get_constraint_type(j);
}

inline ConstraintType ElasticReformulation::get_constraint_status(size_t j) const {
   return this->original_problem.get_constraint_status(j);
}

inline size_t ElasticReformulation::get_hessian_maximum_number_nonzeros() const {
   return this->original_problem.get_hessian_maximum_number_nonzeros();
}

inline void ElasticReformulation::get_initial_primal_point(std::vector<double>& x) const {
   this->original_problem.get_initial_primal_point(x);
   // add the contribution of the elastics
   this->elastic_variables.positive.for_each_value([&](size_t elastic_index) {
      x[elastic_index] = 0.;
   });
   this->elastic_variables.negative.for_each_value([&](size_t elastic_index) {
      x[elastic_index] = 0.;
   });
}

inline void ElasticReformulation::get_initial_dual_point(std::vector<double>& multipliers) const {
   this->original_problem.get_initial_dual_point(multipliers);
}

inline void ElasticReformulation::set_objective_multiplier(double new_objective_multiplier) {
   // update the objective multiplier
   this->objective_multiplier = new_objective_multiplier;
}