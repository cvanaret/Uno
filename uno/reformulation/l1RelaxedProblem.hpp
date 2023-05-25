// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_L1RELAXEDPROBLEM_H
#define UNO_L1RELAXEDPROBLEM_H

#include "RelaxedProblem.hpp"
#include "tools/Range.hpp"
#include "tools/Infinity.hpp"
#include "linear_algebra/VectorExpression.hpp"

struct ElasticVariables {
   SparseVector<size_t> positive;
   SparseVector<size_t> negative;
   explicit ElasticVariables(size_t capacity): positive(capacity), negative(capacity) {}
   [[nodiscard]] size_t size() const { return this->positive.size() + this->negative.size(); }
};

class l1RelaxedProblem: public RelaxedProblem {
public:
   l1RelaxedProblem(const Model& model, double objective_multiplier, double constraint_violation_coefficient);

   [[nodiscard]] double get_objective_multiplier() const override;
   [[nodiscard]] double evaluate_objective(Iterate& iterate) const override;
   void evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const override;
   void evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const override;
   void evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const override;
   void evaluate_lagrangian_hessian(const std::vector<double>& x, const std::vector<double>& multipliers, SymmetricMatrix<double>& hessian) const override;

   void set_infeasibility_measure(Iterate& iterate, Norm progress_norm) const override;
   void set_optimality_measure(Iterate& iterate) const override;
   [[nodiscard]] double compute_predicted_infeasibility_reduction_model(const Iterate& current_iterate, const Direction& direction,
         double step_length, Norm progress_norm) const override;
   [[nodiscard]] std::function<double(double)> compute_predicted_optimality_reduction_model(const Iterate& current_iterate,
         const Direction& direction, double step_length, const SymmetricMatrix<double>& hessian) const override;

   [[nodiscard]] double compute_stationarity_error(const Iterate& iterate, Norm residual_norm) const override;
   [[nodiscard]] double compute_complementarity_error(const std::vector<double>& primals, const std::vector<double>& constraints,
         const Multipliers& multipliers, Norm residual_norm) const override;

   [[nodiscard]] double get_variable_lower_bound(size_t i) const override;
   [[nodiscard]] double get_variable_upper_bound(size_t i) const override;
   [[nodiscard]] double get_constraint_lower_bound(size_t j) const override;
   [[nodiscard]] double get_constraint_upper_bound(size_t j) const override;

   [[nodiscard]] size_t get_number_objective_gradient_nonzeros() const override;
   [[nodiscard]] size_t get_number_jacobian_nonzeros() const override;
   [[nodiscard]] size_t get_number_hessian_nonzeros() const override;

   // parameterization
   void set_objective_multiplier(double new_objective_multiplier);

   void set_elastic_variable_values(Iterate& iterate, const std::function<void(Iterate&, size_t, size_t, double)>& elastic_setting_function) const;

protected:
   double objective_multiplier;
   const double constraint_violation_coefficient;
   ElasticVariables elastic_variables;

   [[nodiscard]] static size_t count_elastic_variables(const Model& model);
   void generate_elastic_variables();
};

inline l1RelaxedProblem::l1RelaxedProblem(const Model& model, double objective_multiplier, double constraint_violation_coefficient):
      RelaxedProblem(model, model.number_variables + l1RelaxedProblem::count_elastic_variables(model), model.number_constraints),
      objective_multiplier(objective_multiplier),
      constraint_violation_coefficient(constraint_violation_coefficient),
      elastic_variables(this->number_constraints) {
   this->generate_elastic_variables();

   // figure out bounded variables
   for (size_t i: this->model.lower_bounded_variables) {
      this->lower_bounded_variables.push_back(i);
   }
   for (size_t i: this->model.upper_bounded_variables) {
      this->upper_bounded_variables.push_back(i);
   }
   for (size_t i: this->model.single_lower_bounded_variables) {
      this->single_lower_bounded_variables.push_back(i);
   }
   for (size_t i: this->model.single_upper_bounded_variables) {
      this->single_upper_bounded_variables.push_back(i);
   }
   const auto add_nonnegative_elastic = [&](size_t elastic_index) {
      this->lower_bounded_variables.push_back(elastic_index);
      this->single_lower_bounded_variables.push_back(elastic_index);
   };
   this->elastic_variables.positive.for_each_value(add_nonnegative_elastic);
   this->elastic_variables.negative.for_each_value(add_nonnegative_elastic);
}

inline double l1RelaxedProblem::get_objective_multiplier() const {
   return this->objective_multiplier;
}

// return rho*f(x) + coeff*||c(x)||₁
inline double l1RelaxedProblem::evaluate_objective(Iterate& iterate) const {
   double objective = 0.;

   // scaled objective: rho*f(x)
   if (this->objective_multiplier != 0.) {
      iterate.evaluate_objective(this->model);
      objective += this->objective_multiplier * iterate.evaluations.objective;
   }

   // scaled constraint violation: coeff*||c(x)||₁
   iterate.evaluate_constraints(this->model);
   objective += this->constraint_violation_coefficient * this->model.compute_constraint_violation(iterate.evaluations.constraints, Norm::L1);
   return objective;
}

inline void l1RelaxedProblem::evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const {
   // scale nabla f(x) by rho
   if (this->objective_multiplier != 0.) {
      iterate.evaluate_objective_gradient(this->model);
      objective_gradient = iterate.evaluations.objective_gradient;
      scale(objective_gradient, this->objective_multiplier);
   }
   else {
      objective_gradient.clear();
   }

   // elastic contribution
   const auto insert_elastic_derivative = [&](size_t elastic_index) {
      objective_gradient.insert(elastic_index, this->constraint_violation_coefficient);
   };
   this->elastic_variables.positive.for_each_value(insert_elastic_derivative);
   this->elastic_variables.negative.for_each_value(insert_elastic_derivative);
}

inline void l1RelaxedProblem::evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const {
   iterate.evaluate_constraints(this->model);
   copy_from(constraints, iterate.evaluations.constraints);
   // add the contribution of the elastics
   this->elastic_variables.positive.for_each([&](size_t j, size_t elastic_index) {
      constraints[j] -= iterate.primals[elastic_index];
   });
   this->elastic_variables.negative.for_each([&](size_t j, size_t elastic_index) {
      constraints[j] += iterate.primals[elastic_index];
   });
}

inline void l1RelaxedProblem::evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const {
   iterate.evaluate_constraint_jacobian(this->model);
   constraint_jacobian = iterate.evaluations.constraint_jacobian;
   // add the contribution of the elastics
   this->elastic_variables.positive.for_each([&](size_t j, size_t elastic_index) {
      constraint_jacobian[j].insert(elastic_index, -1.);
   });
   this->elastic_variables.negative.for_each([&](size_t j, size_t elastic_index) {
      constraint_jacobian[j].insert(elastic_index, 1.);
   });
}

inline void l1RelaxedProblem::evaluate_lagrangian_hessian(const std::vector<double>& x, const std::vector<double>& multipliers,
      SymmetricMatrix<double>& hessian) const {
   this->model.evaluate_lagrangian_hessian(x, this->objective_multiplier, multipliers, hessian);

   // extend the dimension of the Hessian by finalizing the remaining columns (note: the elastics do not enter the Hessian)
   for (size_t j: Range(this->model.number_variables, this->number_variables)) {
      hessian.finalize_column(j);
   }
}

inline void l1RelaxedProblem::set_infeasibility_measure(Iterate& iterate, Norm /*progress_norm*/) const {
   if (this->objective_multiplier == 0.) {
      iterate.progress.infeasibility = 0.;
   }
   else { // 0. < objective_multiplier
      iterate.evaluate_constraints(this->model);
      iterate.progress.infeasibility = this->model.compute_constraint_violation(iterate.evaluations.constraints, Norm::L1);
   }
}

inline void l1RelaxedProblem::set_optimality_measure(Iterate& iterate) const {
   if (this->objective_multiplier == 0.) {
      // constraint violation
      iterate.evaluate_constraints(this->model);
      const double constraint_violation = this->constraint_violation_coefficient *
                                          this->model.compute_constraint_violation(iterate.evaluations.constraints, Norm::L1);
      iterate.progress.optimality = [=](double /*objective_multiplier*/) {
         return constraint_violation;
      };
   }
   else { // 0. < objective_multiplier
      // scaled objective
      iterate.evaluate_objective(this->model);
      const double objective = iterate.evaluations.objective;
      iterate.progress.optimality = [=](double objective_multiplier) {
         return objective_multiplier*objective;
      };
   }
}

// predicted infeasibility reduction
inline double l1RelaxedProblem::compute_predicted_infeasibility_reduction_model(const Iterate& current_iterate, const Direction& direction,
      double step_length, Norm /*progress_norm*/) const {
   if (this->objective_multiplier == 0.) {
      return 0.;
   }
   else { // 0. < objective_multiplier
      // "‖c(x)‖₁ - ‖c(x) + ∇c(x)^T (αd)‖₁"
      const double current_constraint_violation = this->model.compute_constraint_violation(current_iterate.evaluations.constraints, Norm::L1);
      const double trial_linearized_constraint_violation = this->model.compute_linearized_constraint_violation(direction.primals,
            current_iterate.evaluations.constraints, current_iterate.evaluations.constraint_jacobian, step_length, Norm::L1);
      return current_constraint_violation - trial_linearized_constraint_violation;
   }
}

inline std::function<double(double)> l1RelaxedProblem::compute_predicted_optimality_reduction_model(const Iterate& current_iterate,
      const Direction& direction, double step_length, const SymmetricMatrix<double>& hessian) const {
   const double quadratic_problem = hessian.quadratic_product(direction.primals, direction.primals);
   if (this->objective_multiplier == 0.) {
      // "‖c(x)‖₁ - ‖c(x) + ∇c(x)^T (αd)‖₁"
      const double current_constraint_violation = this->model.compute_constraint_violation(current_iterate.evaluations.constraints, Norm::L1);
      const double trial_linearized_constraint_violation = this->model.compute_linearized_constraint_violation(direction.primals,
            current_iterate.evaluations.constraints, current_iterate.evaluations.constraint_jacobian, step_length, Norm::L1);
      return [=](double /*objective_multiplier*/) {
         return this->constraint_violation_coefficient * (current_constraint_violation - trial_linearized_constraint_violation) -
            step_length*step_length/2. * quadratic_problem;
      };
   }
   else { // 0. < objective_multiplier
      // "-ρ*∇f(x)^T (αd)"
      const double directional_derivative = dot(direction.primals, current_iterate.evaluations.objective_gradient);
      return [=](double objective_multiplier) {
         return step_length * (-objective_multiplier*directional_derivative) - step_length*step_length/2. * quadratic_problem;
      };
   }
}

inline double l1RelaxedProblem::compute_stationarity_error(const Iterate& iterate, Norm residual_norm) const {
   // norm of the constraints' contribution of the Lagrangian gradient
   return norm(residual_norm, iterate.lagrangian_gradient.constraints_contribution);
}

// complementary slackness error
inline double l1RelaxedProblem::compute_complementarity_error(const std::vector<double>& primals, const std::vector<double>& constraints,
      const Multipliers& multipliers, Norm residual_norm) const {
   // construct a lazy expression for complementarity for variable bounds
   VectorExpression<double> variable_complementarity(this->get_number_original_variables(), [&](size_t i) {
      if (0. < multipliers.lower_bounds[i]) {
         return multipliers.lower_bounds[i] * (primals[i] - this->model.get_variable_lower_bound(i));
      }
      if (multipliers.upper_bounds[i] < 0.) {
         return multipliers.upper_bounds[i] * (primals[i] - this->model.get_variable_upper_bound(i));
      }
      return 0.;
   });

   // construct a lazy expression for complementarity for constraint bounds
   VectorExpression<double> constraint_complementarity(constraints.size(), [&](size_t j) {
      // violated constraints
      if (constraints[j] < this->get_constraint_lower_bound(j)) { // lower violated
         return (this->constraint_violation_coefficient - multipliers.constraints[j]) * (constraints[j] - this->get_constraint_lower_bound(j));
      }
      else if (this->get_constraint_upper_bound(j) < constraints[j]) { // upper violated
         return (this->constraint_violation_coefficient + multipliers.constraints[j]) * (constraints[j] - this->get_constraint_upper_bound(j));
      }
      // satisfied constraints
      else if (0. < multipliers.constraints[j]) { // lower bound
         return multipliers.constraints[j] * (constraints[j] - this->get_constraint_lower_bound(j));
      }
      else if (multipliers.constraints[j] < 0.) { // upper bound
         return multipliers.constraints[j] * (constraints[j] - this->get_constraint_upper_bound(j));
      }
      return 0.;
   });
   return norm(residual_norm, variable_complementarity, constraint_complementarity);
}

inline double l1RelaxedProblem::get_variable_lower_bound(size_t i) const {
   if (i < this->model.number_variables) { // original variable
      return this->model.get_variable_lower_bound(i);
   }
   else { // elastic variable in [0, +inf[
      return 0.;
   }
}

inline double l1RelaxedProblem::get_variable_upper_bound(size_t i) const {
   if (i < this->model.number_variables) { // original variable
      return this->model.get_variable_upper_bound(i);
   }
   else { // elastic variable in [0, +inf[
      return INF<double>;
   }
}

inline double l1RelaxedProblem::get_constraint_lower_bound(size_t j) const {
   return this->model.get_constraint_lower_bound(j);
}

inline double l1RelaxedProblem::get_constraint_upper_bound(size_t j) const {
   return this->model.get_constraint_upper_bound(j);
}

inline size_t l1RelaxedProblem::get_number_objective_gradient_nonzeros() const {
   size_t number_nonzeros = 0;

   // objective contribution
   if (this->objective_multiplier != 0.) {
      number_nonzeros += this->model.get_number_objective_gradient_nonzeros();
   }

   // elastic contribution
   number_nonzeros += this->elastic_variables.size();
   return number_nonzeros;
}

inline size_t l1RelaxedProblem::get_number_jacobian_nonzeros() const {
   return this->model.get_number_jacobian_nonzeros() + this->elastic_variables.size();
}

inline size_t l1RelaxedProblem::get_number_hessian_nonzeros() const {
   return this->model.get_number_hessian_nonzeros();
}

inline void l1RelaxedProblem::set_objective_multiplier(double new_objective_multiplier) {
   assert(0. <= new_objective_multiplier && "The objective multiplier should be non-negative");
   // update the objective multiplier
   this->objective_multiplier = new_objective_multiplier;
}

inline size_t l1RelaxedProblem::count_elastic_variables(const Model& model) {
   size_t number_elastic_variables = 0;
   // if the subproblem uses slack variables, the bounds of the constraints are [0, 0]
   for (size_t j: Range(model.number_constraints)) {
      if (is_finite(model.get_constraint_lower_bound(j))) {
         number_elastic_variables++;
      }
      if (is_finite(model.get_constraint_upper_bound(j))) {
         number_elastic_variables++;
      }
   }
   return number_elastic_variables;
}

inline void l1RelaxedProblem::generate_elastic_variables() {
   // generate elastic variables to relax the constraints
   size_t elastic_index = this->model.number_variables;
   for (size_t j: Range(this->model.number_constraints)) {
      if (is_finite(this->model.get_constraint_upper_bound(j))) {
         // nonnegative variable p that captures the positive part of the constraint violation
         this->elastic_variables.positive.insert(j, elastic_index);
         elastic_index++;
      }
      if (is_finite(this->model.get_constraint_lower_bound(j))) {
         // nonpositive variable n that captures the negative part of the constraint violation
         this->elastic_variables.negative.insert(j, elastic_index);
         elastic_index++;
      }
   }
}

inline void l1RelaxedProblem::set_elastic_variable_values(Iterate& iterate, const std::function<void(Iterate&, size_t, size_t, double)>&
      elastic_setting_function) const {
   iterate.set_number_variables(this->number_variables);
   this->elastic_variables.positive.for_each([&](size_t j, size_t elastic_index) {
      elastic_setting_function(iterate, j, elastic_index, -1.);
   });
   this->elastic_variables.negative.for_each([&](size_t j, size_t elastic_index) {
      elastic_setting_function(iterate, j, elastic_index, 1.);
   });
}

#endif // UNO_L1RELAXEDPROBLEM_H
