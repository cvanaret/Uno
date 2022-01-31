#ifndef UNO_L1ELASTICREFORMULATION_H
#define UNO_L1ELASTICREFORMULATION_H

#include <vector>
#include <cmath>
#include "Problem.hpp"
#include "Constraint.hpp"
#include "ingredients/constraint_relaxation/ElasticVariables.hpp"

class l1ElasticReformulation: public Problem {
public:
   l1ElasticReformulation(const Problem& original_problem, double objective_multiplier, double elastic_objective_coefficient, bool use_proximal_term);

   [[nodiscard]] size_t get_number_original_variables() const override;
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
   [[nodiscard]] double compute_linearized_constraint_violation(const std::vector<double>& x) const;
   [[nodiscard]] double compute_elastic_residual(const std::vector<double>& x, const std::vector<double>& dx) const;
   [[nodiscard]] double compute_constraint_violation(double constraint, size_t j) const override;

   [[nodiscard]] ConstraintType get_variable_status(size_t i) const override;
   [[nodiscard]] FunctionType get_constraint_type(size_t j) const override;
   [[nodiscard]] ConstraintType get_constraint_status(size_t j) const override;
   [[nodiscard]] size_t get_hessian_maximum_number_nonzeros() const override;

   void get_initial_primal_point(std::vector<double>& x) const override;
   void get_initial_dual_point(std::vector<double>& multipliers) const override;
   void set_objective_multiplier(double new_objective_multiplier);
   void set_proximal_coefficient(double new_proximal_coefficient);
   void set_proximal_reference_point(const std::vector<double>& new_proximal_reference_point);
   void set_elastic_variables(Iterate& first_iterate) const;
   void reset_elastic_variables(Iterate& first_iterate) const;

protected:
   const Problem& original_problem;
   double objective_multiplier;
   // elastic variables
   ElasticVariables elastic_variables;
   double elastic_objective_coefficient;
   // proximal term
   const bool use_proximal_term;
   double proximal_coefficient{0.};
   std::vector<double> proximal_reference_point;

   [[nodiscard]] static size_t count_elastic_variables(const Problem& problem);
   static void generate_elastic_variables(const Problem& problem, ElasticVariables& elastic_variables, size_t number_variables);
   [[nodiscard]] double get_proximal_weight(size_t i) const;
};

inline l1ElasticReformulation::l1ElasticReformulation(const Problem& original_problem, double objective_multiplier,
      double elastic_objective_coefficient, bool use_proximal_term):
      Problem(original_problem.name + "_slacks", // name
            original_problem.number_variables + l1ElasticReformulation::count_elastic_variables(original_problem), // number of variables
            original_problem.number_constraints, // number of constraints
            original_problem.problem_type), // problem type
      original_problem(original_problem),
      objective_multiplier(objective_multiplier),
      // elastic variables
      elastic_variables(this->number_constraints),
      elastic_objective_coefficient(elastic_objective_coefficient),
      use_proximal_term(use_proximal_term),
      proximal_reference_point(original_problem.number_variables) {
   // register equality and inequality constraints
   this->original_problem.equality_constraints.for_each([&](size_t j, size_t i) {
      this->equality_constraints.insert(j, i);
   });
   this->original_problem.inequality_constraints.for_each([&](size_t j, size_t i) {
      this->inequality_constraints.insert(j, i);
   });

   // generate elastic variables
   l1ElasticReformulation::generate_elastic_variables(this->original_problem, this->elastic_variables, this->original_problem.number_variables);

   // figure out bounded variables
   for (size_t i: this->original_problem.lower_bounded_variables) {
      this->lower_bounded_variables.push_back(i);
   }
   for (size_t i: this->original_problem.upper_bounded_variables) {
      this->upper_bounded_variables.push_back(i);
   }
   this->elastic_variables.positive.for_each_value([&](size_t elastic_index) {
      this->lower_bounded_variables.push_back(elastic_index);
   });
   this->elastic_variables.negative.for_each_value([&](size_t elastic_index) {
      this->lower_bounded_variables.push_back(elastic_index);
   });
}

inline size_t l1ElasticReformulation::get_number_original_variables() const {
   return this->original_problem.get_number_original_variables();
}

inline double l1ElasticReformulation::get_variable_lower_bound(size_t i) const {
   if (i < this->original_problem.number_variables) { // original variable
      return this->original_problem.get_variable_lower_bound(i);
   }
   else { // elastic variable in [0, +inf[
      return 0.;
   }
}

inline double l1ElasticReformulation::get_variable_upper_bound(size_t i) const {
   if (i < this->original_problem.number_variables) { // original variable
      return this->original_problem.get_variable_upper_bound(i);
   }
   else { // elastic variable in [0, +inf[
      return std::numeric_limits<double>::infinity();
   }
}

inline double l1ElasticReformulation::get_constraint_lower_bound(size_t j) const {
   return this->original_problem.get_constraint_lower_bound(j);
}

inline double l1ElasticReformulation::get_constraint_upper_bound(size_t j) const {
   return this->original_problem.get_constraint_upper_bound(j);
}

inline double l1ElasticReformulation::compute_linearized_constraint_violation(const std::vector<double>& x) const {
   double residual = 0.;
   // l1 residual of the linearized constraints: sum of elastic variables
   auto elastic_contribution = [&](size_t i) {
      residual += x[i];
   };
   this->elastic_variables.positive.for_each_value(elastic_contribution);
   this->elastic_variables.negative.for_each_value(elastic_contribution);
   assert(0 <= residual && "The elastic residual should not be negative");
   return residual;
}

inline double l1ElasticReformulation::compute_elastic_residual(const std::vector<double>& x, const std::vector<double>& dx) const {
   double residual = 0.;
   // l1 residual of the linearized constraints: sum of elastic variables
   auto elastic_contribution = [&](size_t i) {
      residual += (x[i] + dx[i]);
   };
   this->elastic_variables.positive.for_each_value(elastic_contribution);
   this->elastic_variables.negative.for_each_value(elastic_contribution);
   assert(0 <= residual && "The elastic residual should not be negative");
   return residual;
}

inline double l1ElasticReformulation::compute_constraint_violation(double constraint, size_t j) const {
   return this->original_problem.compute_constraint_violation(constraint, j);
}

// return rho*f(x) + coeff*(e^T p + e^T n) + proximal
inline double l1ElasticReformulation::evaluate_objective(const std::vector<double>& x) const {
   // elastic contribution
   double objective = this->elastic_objective_coefficient* this->compute_linearized_constraint_violation(x);
   // original objective
   if (this->objective_multiplier != 0.) {
      objective += this->objective_multiplier*this->original_problem.evaluate_objective(x);
   }
   // proximal term
   if (this->use_proximal_term && 0. < this->proximal_coefficient) {
      double proximal_term = 0.;
      for (size_t i = 0; i < this->original_problem.number_variables; i++) {
         const double weight = this->get_proximal_weight(i);
         // weighted distance between trial iterate and current iterate
         proximal_term += std::pow(weight * (x[i] - this->proximal_reference_point[i]), 2);
      }
      objective += this->proximal_coefficient*proximal_term;
   }
   return objective;
}

inline void l1ElasticReformulation::evaluate_objective_gradient(const std::vector<double>& x, SparseVector<double>& gradient) const {
   // scale nabla f(x) by rho
   if (this->objective_multiplier != 0.) {
      this->original_problem.evaluate_objective_gradient(x, gradient);
      scale(gradient, this->objective_multiplier);
   }
   else {
      gradient.clear();
   }
   // elastic contribution
   const auto insert_elastic_derivative = [&](size_t elastic_index) {
      gradient.insert(elastic_index, this->elastic_objective_coefficient);
   };
   this->elastic_variables.positive.for_each_value(insert_elastic_derivative);
   this->elastic_variables.negative.for_each_value(insert_elastic_derivative);
   // proximal term
   if (this->use_proximal_term && 0. < this->proximal_coefficient) {
      for (size_t i = 0; i < this->original_problem.number_variables; i++) {
         const double weight = this->get_proximal_weight(i);
         // measure weighted distance between trial iterate and current iterate
         const double derivative = this->proximal_coefficient * weight * (x[i] - this->proximal_reference_point[i]);
         gradient.insert(i, derivative);
      }
   }
}

inline void l1ElasticReformulation::evaluate_constraints(const std::vector<double>& x, std::vector<double>& constraints) const {
   this->original_problem.evaluate_constraints(x, constraints);
   // add the contribution of the elastics
   this->elastic_variables.positive.for_each([&](size_t j, size_t elastic_index) {
      constraints[j] -= x[elastic_index];
   });
   this->elastic_variables.negative.for_each([&](size_t j, size_t elastic_index) {
      constraints[j] += x[elastic_index];
   });
}

inline void l1ElasticReformulation::evaluate_constraint_jacobian(const std::vector<double>& x, std::vector<SparseVector<double>>& constraint_jacobian) const {
   this->original_problem.evaluate_constraint_jacobian(x, constraint_jacobian);
   // add the contribution of the elastics
   this->elastic_variables.positive.for_each([&](size_t j, size_t elastic_index) {
      constraint_jacobian[j].insert(elastic_index, -1.);
   });
   this->elastic_variables.negative.for_each([&](size_t j, size_t elastic_index) {
      constraint_jacobian[j].insert(elastic_index, 1.);
   });
}

inline void l1ElasticReformulation::evaluate_lagrangian_hessian(const std::vector<double>& x, double /*objective_multiplier*/,
      const std::vector<double>& multipliers, SymmetricMatrix& hessian) const {
   this->original_problem.evaluate_lagrangian_hessian(x, this->objective_multiplier, multipliers, hessian);
   hessian.dimension = this->number_variables;
   // add proximal term for the original variables
   if (this->use_proximal_term && 0. < this->proximal_coefficient) {
      for (size_t i = 0; i < this->original_problem.number_variables; i++) {
         const double distance = std::pow(this->get_proximal_weight(i), 2);
         const double diagonal_term = this->proximal_coefficient * distance;
         hessian.insert(diagonal_term, i, i);
      }
   }
   // extend the dimension of the Hessian by finalizing the remaining columns (note: the elastics do not enter the Hessian)
   for (size_t j = this->original_problem.number_variables; j < this->number_variables; j++) {
      hessian.finalize(j);
   }
}

inline ConstraintType l1ElasticReformulation::get_variable_status(size_t i) const {
   if (i < this->original_problem.number_variables) { // original variable
      return this->original_problem.get_variable_status(i);
   }
   else { // elastic variable in [0, +inf[
      return BOUNDED_LOWER;
   }
}

inline FunctionType l1ElasticReformulation::get_constraint_type(size_t j) const {
   return this->original_problem.get_constraint_type(j);
}

inline ConstraintType l1ElasticReformulation::get_constraint_status(size_t j) const {
   return this->original_problem.get_constraint_status(j);
}

inline size_t l1ElasticReformulation::get_hessian_maximum_number_nonzeros() const {
   // add the proximal term
   return this->original_problem.get_hessian_maximum_number_nonzeros() + (use_proximal_term ? this->original_problem.number_variables : 0);
}

inline void l1ElasticReformulation::get_initial_primal_point(std::vector<double>& x) const {
   this->original_problem.get_initial_primal_point(x);
   // add the contribution of the elastics
   this->elastic_variables.positive.for_each_value([&](size_t elastic_index) {
      x[elastic_index] = 0.;
   });
   this->elastic_variables.negative.for_each_value([&](size_t elastic_index) {
      x[elastic_index] = 0.;
   });
}

inline void l1ElasticReformulation::get_initial_dual_point(std::vector<double>& multipliers) const {
   this->original_problem.get_initial_dual_point(multipliers);
}

inline void l1ElasticReformulation::set_objective_multiplier(double new_objective_multiplier) {
   // update the objective multiplier
   this->objective_multiplier = new_objective_multiplier;
}

inline void l1ElasticReformulation::set_proximal_coefficient(double new_proximal_coefficient) {
   // update the proximal coefficient
   this->proximal_coefficient = new_proximal_coefficient;
}

inline void l1ElasticReformulation::set_proximal_reference_point(const std::vector<double>& new_proximal_reference_point) {
   // update the proximal reference point
   copy_from(this->proximal_reference_point, new_proximal_reference_point, this->original_problem.number_variables);
}

inline size_t l1ElasticReformulation::count_elastic_variables(const Problem& problem) {
   size_t number_elastic_variables = 0;
   // if the subproblem uses slack variables, the bounds of the constraints are [0, 0]
   for (size_t j = 0; j < problem.number_constraints; j++) {
      if (is_finite(problem.get_constraint_lower_bound(j))) {
         number_elastic_variables++;
      }
      if (is_finite(problem.get_constraint_upper_bound(j))) {
         number_elastic_variables++;
      }
   }
   return number_elastic_variables;
}

inline void l1ElasticReformulation::generate_elastic_variables(const Problem& problem, ElasticVariables& elastic_variables, size_t number_variables) {
   // generate elastic variables p and n on the fly to relax the constraints
   // if the subproblem uses slack variables, the bounds of the constraints are [0, 0]
   size_t elastic_index = number_variables;
   for (size_t j = 0; j < problem.number_constraints; j++) {
      if (is_finite(problem.get_constraint_upper_bound(j))) {
         // nonnegative variable p that captures the positive part of the constraint violation
         elastic_variables.positive.insert(j, elastic_index);
         elastic_index++;
      }
      if (is_finite(problem.get_constraint_lower_bound(j))) {
         // nonpositive variable n that captures the negative part of the constraint violation
         elastic_variables.negative.insert(j, elastic_index);
         elastic_index++;
      }
   }
}

inline double l1ElasticReformulation::get_proximal_weight(size_t i) const {
   // weight of each diagonal term of the proximal term
   return std::min(1., 1. / std::abs(this->proximal_reference_point[i]));
}

inline void l1ElasticReformulation::set_elastic_variables(Iterate& first_iterate) const {
   first_iterate.set_number_variables(this->number_variables);
   this->elastic_variables.positive.for_each([&](size_t j, size_t elastic_index) {
      first_iterate.x[elastic_index] = this->compute_constraint_upper_bound_violation(first_iterate.constraints[j], j);
   });
   this->elastic_variables.negative.for_each([&](size_t j, size_t elastic_index) {
      first_iterate.x[elastic_index] = this->compute_constraint_lower_bound_violation(first_iterate.constraints[j], j);
   });
}

inline void l1ElasticReformulation::reset_elastic_variables(Iterate& first_iterate) const {
   first_iterate.set_number_variables(this->number_variables);
   this->elastic_variables.positive.for_each_value([&](size_t elastic_index) {
      first_iterate.x[elastic_index] = 0.;
   });
   this->elastic_variables.negative.for_each_value([&](size_t elastic_index) {
      first_iterate.x[elastic_index] = 0.;
   });
}

#endif // UNO_L1ELASTICREFORMULATION_H